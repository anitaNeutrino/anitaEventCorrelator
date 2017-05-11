////////////////////////////////////////
//  BedmapReader.cxx : 
//
//  Based on code stolen from Stephen Hoover's readBedmap.C.  
//  Modified by Matthew Mottram to work for surface elevation
//  and to be consistent with the eventCorrelator framework
//
//    Sample code for reading in BEDMAP data and making it available 
//  by longitude and latitude coordinates.  BEDMAP data is assumed to be 
//  in directory "data/".  The files "icethic.asc", "groundbed.asc", and 
//  "water.asc" are needed.
//    Methods "IceThickness", "Surface", "SurfaceAboveGeoid", and "WaterDepth"
//  take (latitude, longitude) as arguments, and return data from BEDMAP.
//  Longitude is in degrees from the Prime Meridian, with degrees east being 
//  positive and degrees west being negative.  Latitude is in degrees from 
//  the equator, where negative degrees are south latitude.
//
//    Code to read data from files by Ryan Nichol, remaining code by Stephen Hoover.
//  16 April 2007
//
//////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include "TMath.h"
#include "BedmapReader.h"
#include "AnitaGeomTool.h"
#include "TProfile2D.h"
#include "TGaxis.h"

using namespace std;

 struct rampData{
  Double_t lat;
  Double_t longi;
  Double_t el1;
  Double_t el2;
 };

bool comparitor(rampData const& first, rampData const& second){
	return first.lat < second.lat;
}

vector<rampData> rampDemVec;
vector<double> longVec;
vector<double> latVec;

vector<double>::iterator myIteratorLow;
vector<double>::iterator myIteratorUp;
vector<double>::iterator myIteratorFinal;

//Parameters of the BEDMAP ice model. (See http://www.antarctica.ac.uk/data/access/bedmap/download/)
// Int_t cellSize=5000; //in meters, set by header file (should be 5000)
// Int_t nCols_surface=1096;
// Int_t nRows_surface=911;
// Int_t xLowerLeft_surface=-2713600;
// Int_t yLowerLeft_surface=-2304000;
Int_t NODATA=-9999;

//Variables for conversion between BEDMAP polar stereographic coordinates and lat/lon.  Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf
const Double_t scale_factor=0.97276901289;  //scale factor at pole corresponding to 71 deg S latitude of true scale (used in BEDMAP)
const Double_t ellipsoid_inv_f = 298.257223563; //of Earth
const Double_t ellipsoid_b = R_EARTH*(1-(1/ellipsoid_inv_f));
const Double_t eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));
const Double_t bedmap_a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;
const Double_t bedmap_b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;
const Double_t bedmap_c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;
const Double_t bedmap_d_bar = 4279*pow(eccentricity,8)/161280;
const Double_t bedmap_c_0 = (2*R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);
Double_t bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(71*TMath::DegToRad())) / (1 - eccentricity*sin(71*TMath::DegToRad())) ),eccentricity/2) * tan((TMath::Pi()/4) - (71*TMath::DegToRad())/2); //varies with latitude, defined here for 71 deg S latitude
const Double_t bedmap_nu = bedmap_R / cos(71*TMath::DegToRad());


BedmapReader*  BedmapReader::fgInstance = 0;

BedmapReader::BedmapReader(bool icethicknessMode) 
{
  //Default constructor
  std::cout << "reading the bedmap data" << std::endl;
  ReadSurfaceElevation(icethicknessMode);
  //ReadSurfaceElevationRampDem();
  fgInstance=this;
}


BedmapReader::~BedmapReader() 
{
  //Default destructor
}


//______________________________________________________________________________
BedmapReader*  BedmapReader::Instance(bool icethicknessMode)
{
  //static function
  return (fgInstance) ? (BedmapReader*) fgInstance : new BedmapReader(icethicknessMode);
}



//_______________________________________________________________________________
Double_t BedmapReader::Surface(Double_t lon,Double_t lat) {
  return (SurfaceAboveGeoid(lon,lat) + Geoid(lat));
} //Surface




//_______________________________________________________________________________
Double_t BedmapReader::SurfaceAboveGeoid(Double_t lon, Double_t lat) {
  //This method returns the elevation above the geoid of the surface of the top of the ice ice (or bare ground, if no ice is present) in meters, at a location specified by a latitude and longitude (in degrees). 
  Double_t surface=0;

  Int_t e_coord_surface=0;
  Int_t n_coord_surface=0;
  SurfaceLonLattoEN(lon,lat,e_coord_surface,n_coord_surface);

  if (e_coord_surface <= nCols_surface && e_coord_surface >= 0 && n_coord_surface <= nRows_surface && n_coord_surface >= 0)
    surface = surface_elevation[e_coord_surface][n_coord_surface];

  else {
    surface = NODATA;
  }
//   std::cout << "surface height " << surface << " e_co " << e_coord_surface << " n_co " << n_coord_surface << " lon " << lon << " lat " << lat << std::endl;
  return surface;
} //method SurfaceAboveGeoid



//_______________________________________________________________________________
Double_t BedmapReader::Geoid(Double_t latitude) {
  return (GEOID_MIN*GEOID_MAX/sqrt(pow(GEOID_MIN,2)-(pow(GEOID_MIN,2)-pow(GEOID_MAX,2))*pow(cos(latitude*TMath::DegToRad()),2)));
} //Geoid(lat)

//_______________________________________________________________________________
Double_t BedmapReader::SurfaceAboveGeoidRampDem(Double_t lon, Double_t lat) {

  double searchVal = lat;
  double searchValFinal = lon;
  double lower = 0;
  double  upper = 0;
  double finalVal = 0;
  double lowerEl =0;
  double upperEl =0;
  double finalEl;


  //there may be more than one value that matches the latitude, so first lets fins find the upper and lower limit on this range
 myIteratorLow = lower_bound(latVec.begin(), latVec.end(), searchVal);
if (myIteratorLow == latVec.begin()) upper = *myIteratorLow; // no smaller value  than val in vector
else if (myIteratorLow == latVec.end()) lower = *(myIteratorLow-1); // no bigger value than val in vector
else {
    lower = *(myIteratorLow);
   lowerEl = distance(latVec.begin(),myIteratorLow); 
}

 myIteratorUp = upper_bound(latVec.begin(), latVec.end(), searchVal);
if (myIteratorUp == latVec.begin()) upper = *myIteratorUp; // no smaller value  than val in vector
else if (myIteratorUp == latVec.end()) lower = *(myIteratorUp-1); // no bigger value than val in vector
else {
    upper = *(myIteratorUp);
   upperEl = distance(latVec.begin(),myIteratorUp); 
}


//now we have upper and lower limit, lets find the nearest longitude in this range of values
 myIteratorFinal = lower_bound(longVec.begin()+lowerEl, longVec.begin()+upperEl-1, searchValFinal);
    
   finalVal = *(myIteratorFinal);
   finalEl = distance(longVec.begin()+lowerEl,myIteratorFinal); 

    finalEl = rampDemVec[finalEl+lowerEl].el1;
    //cout << finalEl << endl;

  return finalEl;


}


//_______________________________________________________________________________
void BedmapReader::ReadSurfaceElevation(bool icethicknessMode) {
  //Reads the BEDMAP data on the elevation of the surface beneath the ice.  If there is water beneath the ice, the ground elevation is given the value 0.  Assumes the file is in directory "data".  Origianl code by Ryan Nichol.
  char calibDir[FILENAME_MAX];
  char *calibEnv=getenv("ANITA_CALIB_DIR");
  if(!calibEnv) {
     char *utilEnv=getenv("ANITA_UTIL_INSTALL_DIR");
     if(!utilEnv)
	sprintf(calibDir,"calib");
     else
	sprintf(calibDir,"%s/share/anitaCalib",utilEnv);
  }
  else {
    strncpy(calibDir,calibEnv,FILENAME_MAX);
  }
  char surfaceFile[FILENAME_MAX];
  if(!icethicknessMode)
    sprintf(surfaceFile,"%s/surfaceElevation.asc",calibDir);
  else 
    sprintf(surfaceFile,"/home/mottram/work/eventCorrelator/data/iceThickness.asc");
  ifstream SurfaceElevationFile(surfaceFile);
  if(!SurfaceElevationFile) {
    std::cerr << "Couldn't open: " << surfaceFile << std::endl;
    exit(1);
  }

  std::cout << calibDir << " " << surfaceFile << std::endl;

  std::cout<<"Reading in BEDMAP data on surface elevation.\n";

  string tempBuf1;
  string tempBuf2;
  string tempBuf3;
  string tempBuf4;
  string tempBuf5;
  string tempBuf6;
  Int_t temp1,temp2,temp3,temp4,temp5,temp6;
  
  SurfaceElevationFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		>> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		>> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
  
  if(tempBuf1 == string("ncols")) {
    nCols_surface=temp1;
  }
  if(tempBuf2 == string("nrows")) {
    nRows_surface=temp2;
  }
  if(tempBuf3 == string("xllcorner")) {
    xLowerLeft_surface=temp3;
  }
  if(tempBuf4 == string("yllcorner")) {
    yLowerLeft_surface=temp4;
  }
  if(tempBuf5 == string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == string("NODATA_value")) {
    NODATA=temp6;
  }

  xUpperRight_surface = xLowerLeft_surface + cellSize*nCols_surface;
  yUpperRight_surface = yLowerLeft_surface + cellSize*nRows_surface;

  //std::cout<<"nCols_surface, nRows_surface "<<nCols_surface<<" , "<<nRows_surface<<std::endl;
  //std::cout<<"xLL_surface, yLL_surface, cellsize "<<xLowerLeft_surface<<" , "<<yLowerLeft_surface<<" , "<<cellSize<<std::endl<<std::endl;
  
  surface_elevation = std::vector< std::vector<short> >( nCols_surface, std::vector<short>( nRows_surface, 0 ) );

  Double_t theValue;
  for(Int_t rowNum=0;rowNum<nRows_surface;rowNum++) {
    for(Int_t colNum=0;colNum<nCols_surface;colNum++) {

      SurfaceElevationFile >> theValue;

      
      // if(theValue==NODATA)
      //	theValue=0; //Set elevation to 0 where we have no data.

//       surface_elevation[colNum][rowNum] = Double_t(theValue);
      surface_elevation[colNum][rowNum] = theValue;

      //if (theValue != -96 && theValue != 0)
      //std::cout<<"surface_elevation: "<<theValue<<std::endl;
    }//for
  }//for
  
  SurfaceElevationFile.close();
  return;
} //method ReadSurfaceElevation


//_______________________________________________________________________________
void BedmapReader::ReadSurfaceElevationRampDem() {
  //Reads the BEDMAP data on the elevation of the surface beneath the ice.  If there is water beneath the ice, the ground elevation is given the value 0.  Assumes the file is in directory "data".  Origianl code by Ryan Nichol.
  char calibDir[FILENAME_MAX];
  char *calibEnv=getenv("ANITA_CALIB_DIR");
  if(!calibEnv) {
     char *utilEnv=getenv("ANITA_UTIL_INSTALL_DIR");
     if(!utilEnv)
	sprintf(calibDir,"calib");
     else
	sprintf(calibDir,"%s/share/anitaCalib",utilEnv);
  }
  else {
    strncpy(calibDir,calibEnv,FILENAME_MAX);
  }
  char surfaceFile[FILENAME_MAX];
  //sprintf(surfaceFile,"%s/ramp1kmdem_wgsosu_v2.txt",calibDir);
  sprintf(surfaceFile,"/unix/anita1/rampDemData/ramp1kmdem_wgsosu_v2.txt");
  ifstream SurfaceElevationFile(surfaceFile);
  if(!SurfaceElevationFile) {
    std::cerr << "Couldn't open: " << surfaceFile << std::endl;
    exit(1);
  }

  std::cout<<"Reading in RampDEM data on surface elevation.\n";

 
  Double_t latIn = 0;
  Double_t longIn = 0;
  Double_t el1In = 0;
  Double_t el2In = 0;  


  rampData dataIn;

  while (!SurfaceElevationFile.eof()){
	SurfaceElevationFile >> latIn;
	SurfaceElevationFile >> longIn;
       	SurfaceElevationFile >> el1In;
       	SurfaceElevationFile >> el2In;

	dataIn.lat = latIn;
	dataIn.longi = longIn;
	dataIn.el1 = el1In;
	dataIn.el2 = el2In;

	rampDemVec.push_back(dataIn);

       }

  cout << "sorting data" << endl;
  std::sort(rampDemVec.begin(), rampDemVec.end(), comparitor);

 for (vector<rampData>::iterator it = rampDemVec.begin(); it!=rampDemVec.end(); ++it) {
 
   latVec.push_back(it->lat);
   longVec.push_back(it->longi);
   //cout << it->lat << "  " << it->longi << endl;

  }

  
  SurfaceElevationFile.close();
  return;
} //method ReadSurfaceElevation





//_______________________________________________________________________________
Double_t BedmapReader::Area(Double_t latitude) {
  //Returns the area of one square of the BEDMAP data at a given latitude. 
  Double_t lat_rad = -latitude * TMath::DegToRad();
  //Double_t lat_rad = (90 - latitude) * TMath::DegToRad();

  return (pow(cellSize* ((1 + sin(71*TMath::DegToRad())) / (1 + sin(lat_rad))),2));
} //method Area




//_______________________________________________________________________________
void BedmapReader::LonLattoEN(Double_t lon, Double_t lat, Double_t xLowerLeft, Double_t yLowerLeft, int& e_coord, int& n_coord) {
  //takes as input a latitude and longitude (in degrees) and converts to indicies for BEDMAP matricies. Needs a location for the corner of the matrix, as not all the BEDMAP files cover the same area.  Code by Stephen Hoover.

  Double_t easting=0;
  Double_t northing=0;

  Double_t lon_rad = lon * TMath::DegToRad(); //convert to radians
  Double_t lat_rad = -lat * TMath::DegToRad();

  bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(lat_rad)) / (1 - eccentricity*sin(lat_rad)) ),eccentricity/2) * tan((TMath::Pi()/4) - lat_rad/2);

  easting = bedmap_R * sin(lon_rad);
  northing = bedmap_R * cos(lon_rad);

  e_coord = (int)((easting - xLowerLeft) / cellSize);
  n_coord = (int)((-1*northing - yLowerLeft) / cellSize);

  return;
} //method LonLattoEN




//_______________________________________________________________________________
void BedmapReader::SurfaceLonLattoEN(Double_t lon, Double_t lat, int& e_coord, int& n_coord) {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP surface elevation data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_surface, yLowerLeft_surface, e_coord, n_coord);
}//SurfaceLonLattoEN



//_______________________________________________________________________________
void BedmapReader::ENtoLonLat(Int_t e_coord, Int_t n_coord, Double_t xLowerLeft, Double_t yLowerLeft, Double_t& lon, Double_t& lat) {
  //Takes as input the indicies from a BEDMAP data set, and turns them into latitude and longitude coordinates.  Information on which data set (surface data, ice depth, water depth) is necessary, in the form of coordinates of a corner of the map.  Code by Stephen Hoover.

  Double_t isometric_lat=0;
  Double_t easting = xLowerLeft+(cellSize*(e_coord+0.5)); //Add offset of 0.5 to get coordinates of middle of cell instead of edges.
  Double_t northing = -1*(yLowerLeft+(cellSize*(n_coord+0.5)));

  //first set longitude

  if (northing!=0)
    lon = atan(easting/northing);
  else
    lon = 90*TMath::DegToRad();


  if (easting > 0 && lon < 0) //adjust sign of longitude
    lon += TMath::Pi();
  else if (easting < 0 && lon > 0)
    lon -= TMath::Pi();
  else if (easting == 0 && northing < 0)
    lon += TMath::Pi();

  //now find latitude

  if (easting != 0)
    bedmap_R = TMath::Abs(easting/sin(lon));
  else if (easting == 0 && northing != 0)
    bedmap_R = TMath::Abs(northing);
  else {
    lat = 0; //at the pole, set lat=0 degrees
    lon = lon*TMath::RadToDeg();
    return;
  } //else

  isometric_lat = (TMath::Pi()/2) - 2*atan(bedmap_R/(scale_factor*bedmap_c_0));

  lat = isometric_lat + bedmap_a_bar*sin(2*isometric_lat) + bedmap_b_bar*sin(4*isometric_lat) + bedmap_c_bar*sin(6*isometric_lat) + bedmap_d_bar*sin(8*isometric_lat);

  lon = lon * TMath::RadToDeg();  //convert to degrees
  lat =  -lat*TMath::RadToDeg(); //convert to degrees, with -90 degrees at the south pole
  return;
} //method ENtoLonLat




//_______________________________________________________________________________
void BedmapReader::SurfaceENtoLonLat(Int_t e, Int_t n, Double_t& lon, Double_t& lat) {
  //Converts indicies of the BEDMAP surface elevation matrix into longitude and latitude.  Code by Stephen Hoover.
  ENtoLonLat(e,n,xLowerLeft_surface,yLowerLeft_surface,lon,lat);
}//SurfaceENtoLonLat




//_______________________________________________________________________________

TProfile2D *BedmapReader::bedmapMap(int coarseness_factor, int set_log_scale,UInt_t &xBins,UInt_t &yBins){

  Bool_t debug=false;

  UInt_t num_columns = surface_elevation.size()/coarseness_factor;
  UInt_t num_rows = surface_elevation[0].size()/coarseness_factor;
  //  UInt_t num_columns = 1096/coarseness_factor;
  //  UInt_t num_rows = 911/coarseness_factor;
  xBins = num_columns;
  yBins = num_rows;

  std::cout << "xBins " << xBins << " yBins " << yBins << std::endl;

//   TProfile2D * theHist = new TProfile2D( "antarctica_surface_elevation_hist", "", num_columns, x_min, x_max, num_rows, y_min, y_max+cell_size );
  TProfile2D * theHist = new TProfile2D( "antarctica_ice_thickness_hist", "", num_columns, xLowerLeft_surface, xUpperRight_surface, num_rows, yLowerLeft_surface, yUpperRight_surface );

  for (unsigned int row_index=0; row_index < surface_elevation[0].size(); ++row_index){
    if(debug){
      if(row_index%(surface_elevation[0].size()/100)==0) std::cerr << "*";
    }
	
    for (unsigned int column_index=0; column_index < surface_elevation.size(); ++column_index)
      {

	if(surface_elevation[column_index][row_index]<-9000)//this is a silly scale for map drawing
	  surface_elevation[column_index][row_index]=-500;

	if (set_log_scale == 1)
	  theHist->Fill( xLowerLeft_surface + double(column_index)*cellSize, -(yLowerLeft_surface + double(row_index)*cellSize), log10(surface_elevation[column_index][row_index]+1000) ); //The "+1000" makes sure everything is positive.
// 	  theHist->Fill( x_min + double(column_index)*cell_size, -(y_min + double(row_index)*cell_size), log10(surface_elevation[column_index][row_index]+1000) ); //The "+1000" makes sure everything is positive.
	else
	  theHist->Fill( xLowerLeft_surface + double(column_index)*cellSize, -(yLowerLeft_surface + double(row_index)*cellSize), surface_elevation[column_index][row_index] );
// 	  theHist->Fill( x_min + double(column_index)*cell_size, -(y_min + double(row_index)*cell_size), surface_elevation[column_index][row_index] );
      } //end for (loop over surface elevation vector)
  }
  if(debug)
    std::cout << std::endl;

  theHist->SetStats(0);

  return theHist;


}




//_______________________________________________________________________________

TGaxis *BedmapReader::distanceScale(Double_t xMin,Double_t xMax,Double_t yMin,Double_t yMax){
  //here xmin, xmax etc are the positions in easting/northing of the distance scale
  if(xMax<=xMin && yMax<=yMin) {
    std::cerr << "size ordering wrong: xMin " << xMin << " xMax " << xMax << " yMin " << yMin << " yMax " << yMax << std::endl;
    return NULL;
  }

  TGaxis *theAxis = new TGaxis(xMin,yMin,xMax,yMax,0,sqrt((xMax-xMin)/1e3*(xMax-xMin)/1e3-(yMax-yMin)/1e3*(yMax-yMin)/1e3),2,"");
  theAxis->SetTitle("km");

  return theAxis;
}
