////////////////////////////////////////
//  RampdemReader.cxx :
//
//  More code stolen from Stephen's Antarctica.cxx
//  that will read in Rampdem data to use with
//  UsefulAdu5Pat.cxx to locate sources on the continent
//
//
//////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include "TMath.h"
#include "RampdemReader.h"
#include "AnitaGeomTool.h"
#include "TProfile2D.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TColor.h"



// our map type
typedef std::vector<std::vector<short> > VecVec;
typedef std::map<RampdemReader::dataSet, VecVec > DataMap;
static DataMap bedMap2Data;
typedef std::map<RampdemReader::dataSet, Int_t> HeaderMap;

static HeaderMap numXs;
static HeaderMap numYs;
static HeaderMap noDatas;
static HeaderMap minXs;
static HeaderMap minYs;
static HeaderMap maxXs;
static HeaderMap maxYs;
static HeaderMap cellSizes;




//Variables for conversion between polar stereographic coordinates and lat/lon.  Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf
static double scale_factor=0.97276901289;  //scale factor at pole corresponding to 71 deg S latitude of true scale (used in both BEDMAP and RAMP DEM)
static double ellipsoid_inv_f = 298.257223563; //of Earth
static double ellipsoid_b = R_EARTH*(1-(1/ellipsoid_inv_f));
static double eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));
static double a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;
static double b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;
static double c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;
static double d_bar = 4279*pow(eccentricity,8)/161280;
static double c_0 = (2*R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);
static double R_factor = scale_factor*c_0 * pow(( (1 + eccentricity*sin(71*TMath::RadToDeg())) / (1 - eccentricity*sin(71*TMath::RadToDeg())) ),eccentricity/2) * tan((TMath::Pi()/4) - (71*TMath::RadToDeg())/2); //varies with latitude, defined here for 71 deg S latitude
static double nu_factor = R_factor / cos(71*TMath::RadToDeg());



RampdemReader*  RampdemReader::fgInstance = 0;




/** RAMP DEM data.  Note: x increases to the right, y increases downward.  **/
static std::vector< std::vector<short> > surface_elevation;
static double cell_size;
static double x_min;
static double x_max;
static double y_min;
static double y_max;
static int nRows_surface;
static int  nCols_surface;
static int  nBytes_surface;




RampdemReader::RampdemReader()
{
  //Default constructor
  // std::cout << "reading the rampdem data" << std::endl;
  readRAMPDEM();
  fgInstance=this;
}


RampdemReader::~RampdemReader()
{
  //Default destructor
}

RampdemReader*  RampdemReader::Instance()
{
  //static function
  return (fgInstance) ? (RampdemReader*) fgInstance : new RampdemReader();
}



//_______________________________________________________________________________
Double_t RampdemReader::Surface(Double_t lon,Double_t lat) {
  return (SurfaceAboveGeoid(lon,lat) + Geoid(lat));
} //Surface



/**
 * Returns the elevation above the geoid of the surface of the top of the ice (or bare ground if no ice)
 * in meters, at a location specified by a latitude and longitude (in degrees).
 * @param lon is the longiutde (degrees)
 * @param lat is the latitude (degrees)
 *
 * @return elevation above geoid in metres
 */
Double_t RampdemReader::SurfaceAboveGeoid(Double_t lon, Double_t lat) {

  Double_t surface=0;

  Int_t e_coord_surface=0;
  Int_t n_coord_surface=0;
  LonLattoEN(lon,lat,e_coord_surface,n_coord_surface);

  if(e_coord_surface >= nCols_surface || e_coord_surface <0){
//     std::cerr<<"[RampdemReader::surfaceAboveGeoid]  Error!  Trying to access x-element "<<e_coord_surface<<" of the RAMP DEM data! (Longitude, latitude = "<<lon<<", "<<lat<<")\n";
    return -9999;
  }
  else if(n_coord_surface >= nRows_surface || n_coord_surface <0){
    //     std::cerr<<"[RampdemReader::surfaceAboveGeoid]  Error!  Trying to access y-element "<<n_coord_surface<<" of the RAMP DEM data! (Longitude, latitude = "<<lon<<", "<<lat<<")\n";
    return -9999;
  }
  else{
    surface = double(surface_elevation[e_coord_surface][n_coord_surface]);
  }

  return surface;
}



//_______________________________________________________________________________
Double_t RampdemReader::Geoid(Double_t latitude) {
  return (GEOID_MIN*GEOID_MAX/sqrt(pow(GEOID_MIN,2)-(pow(GEOID_MIN,2)-pow(GEOID_MAX,2))*pow(cos(latitude*TMath::DegToRad()),2)));
} //Geoid(lat)






//_______________________________________________________________________________
int RampdemReader::readRAMPDEM(){
  bool debug=false;

  char calibDir[FILENAME_MAX];
  char *calibEnv=getenv("ANITA_CALIB_DIR");
  if(!calibEnv) {
    char *utilEnv=getenv("ANITA_UTIL_INSTALL_DIR");
    if(!utilEnv){
      sprintf(calibDir,"calib");
    }
    else{
      sprintf(calibDir,"%s/share/anitaCalib",utilEnv);
    }
  }
  else {
    strncpy(calibDir,calibEnv,FILENAME_MAX);
  }

  char dem_filename[FILENAME_MAX];
  char header_filename[FILENAME_MAX];

  std::ifstream dem_data;
  std::ifstream dem_header;

  sprintf(dem_filename,"%s/ramp1kmdem_wgs_v2.bin",calibDir);
  sprintf(header_filename,"%s/ramp1kmdem_wgs_v2.hdr",calibDir);

  /*
    Open header and binary files.
    Check to make sure the opening was successful, and return
    error code 1 if it wasn't.
  */

  dem_header.open(header_filename);
  if (!dem_header.is_open()){
    std::cerr << "[RampdemReader::readRAMPDEM] Error! Could not open header file " << header_filename
	      << "! Exiting method.\n";
    return 1;
  }

  dem_data.open(dem_filename, std::ios::in | std::ios::binary );
  if (!dem_data.is_open()){
    std::cerr << "[RampdemReader::readRAMPDEM] Error! Could not open binary data file " << dem_filename
	      << "! Exiting method.\n";
    return 1;
  }


  // Read the header.
  // An equals sign preceeds each interesting number, so we can ignore everything up to that.

  // int nRows_surface, nCols_surface, nBytes_surface;
  double min_value, max_value, mean, std_deviation;
  dem_header.ignore( 5000, '=');
  dem_header >> cell_size;
  dem_header.ignore( 5000, '=');
  dem_header >> nRows_surface;
  dem_header.ignore( 5000, '=');
  dem_header.ignore( 5000, '=');
  dem_header >> nCols_surface;
  dem_header.ignore( 5000, '=');
  dem_header >> nBytes_surface;
  dem_header.ignore( 5000, '=');
  dem_header >> x_min;
  dem_header.ignore( 5000, '=');
  dem_header >> min_value;
  dem_header.ignore( 5000, '=');
  dem_header >> x_max;
  dem_header.ignore( 5000, '=');
  dem_header >> max_value;
  dem_header.ignore( 5000, '=');
  dem_header >> y_min;
  dem_header.ignore( 5000, '=');
  dem_header >> mean;
  dem_header.ignore( 5000, '=');
  dem_header >> y_max;
  dem_header.ignore( 5000, '=');
  dem_header >> std_deviation;

  numXs[RampdemReader::rampdem] = nRows_surface;
  numYs[RampdemReader::rampdem] = nCols_surface;
  minXs[RampdemReader::rampdem] = x_min;
  minYs[RampdemReader::rampdem] = y_min;
  maxXs[RampdemReader::rampdem] = x_max;
  maxYs[RampdemReader::rampdem] = y_max;

  if (debug){
    std::cout << "cell size = " << cell_size
	      << ", x_min = " << x_min << ", x_max = " << x_max
	      << ", y_min = " << y_min << ", y_max = " << y_max
	      << ", mean = " << mean << std::endl;
  }

  /* Now that we know the size of the grid, allocate the memory to store it. */
  surface_elevation = std::vector< std::vector<short> >(nCols_surface, std::vector<short>(nRows_surface, 0 ) );

  /* Read in the data and store it to the vector. */
  short temp_data=0;
  for (unsigned int row_index=0; row_index < surface_elevation[0].size(); ++row_index){
    for (unsigned int column_index=0; column_index < surface_elevation.size(); ++column_index){
      if (!dem_data.read( (char *)&temp_data, sizeof(short) )){
	std::cerr << "[RampdemReader::readRAMPDEM] Error! Read from data file failed! Row " << row_index
		  << ", column " << column_index << ".\n";
	return 2;
      } //end if
	// RampdemReader::flipEndian( temp_data );
      flipEndian(temp_data);
      surface_elevation[column_index][row_index] = temp_data;
    } //end for
  }
  if(dem_data.read((char*)&temp_data, sizeof(short))){
    std::cerr << "[RampdemReader::readRAMPDEM] Error! I could read a value (" << temp_data
	      << ") from the data file after I thought that it was empty!.\n";
    return 2;
  }

  // For debugging purposes, look through the data we read in.
  if (debug){
    double my_mean = 0;
    int entries = 0;
    short my_min=1;
    short my_max=1;
    for (unsigned int row_index=0; row_index < surface_elevation[0].size(); ++row_index){
      for (unsigned int column_index=0; column_index < surface_elevation.size(); ++column_index){
	int test_int = surface_elevation[column_index][row_index];
	//if (test_int != 0)
	{
	  my_mean += double(test_int);
	  ++entries;
	} //end if
	if (test_int < my_min) my_min = test_int;
	if (test_int > my_max) my_max = test_int;
      } //end while
      std::cout << entries << " entries, mean is " << my_mean/double(entries) << std::endl;
      std::cout << "I found minimum value " << my_min << " and maximum value " << my_max << std::endl;
    } //end if (debug)
  }

  dem_header.close();
  dem_data.close();

  return 0;
}




/**
 * Returns the area of one square of the BEDMAP data at a given latitude.
 *
 * @param latitude is in degrees
 *
 * @return the area of in metres
 */
Double_t RampdemReader::Area(Double_t latitude) {

  Double_t lat_rad = -latitude * TMath::DegToRad();

  return (pow(cell_size* ((1 + sin(71*TMath::DegToRad())) / (1 + sin(lat_rad))),2));
}



/**
 * Takes latitude and longitude (in degrees) and converts to indicies for BEDMAP matricies.
 * Needs a location for the corner of the matrix, as not all the BEDMAP files cover the same area.
 * Code by Stephen Hoover.
 *
 * @param lon is the longitude in degrees
 * @param lat is the latitude in degrees
 * @param e_coord is the bedmap east coordinate index
 * @param n_coord is the bedmap north coordinate index
 */
void RampdemReader::LonLattoEN(Double_t lon, Double_t lat, int& e_coord, int& n_coord) {
  //
  bool debug=false;

  Double_t easting=0;
  Double_t northing=0;

  LonLatToEastingNorthing(lon,lat,easting,northing);
  EastingNorthingToEN(easting,northing,e_coord,n_coord);

  if(debug){
    std::cout << "lon " << lon << " lat " << lat << " easting " << easting << " northing " << northing << " e_coord " << e_coord << " n_coord " << n_coord << std::endl;
  }
}




/**
 * Converts Easting/northing to RAMPDEM data indices
 *
 * @param easting
 * @param northing
 * @param e_coord is the easting coordinate index
 * @param n_coord is the the northing coordinate index
 */
void RampdemReader::EastingNorthingToEN(Double_t easting,Double_t northing,Int_t &e_coord,Int_t &n_coord){
  bool debug=false;
  if(debug){
    std::cout << "easting " << easting << " northing " << northing << " e_coord " << e_coord << " n_coord " << n_coord << " x_min " << x_min << " y_min " << y_min << " cell size " << cell_size << std::endl;
  }
  e_coord = (int)((easting - x_min) / cell_size);
  n_coord = (int)((-1*northing - y_min) / cell_size);

  if(debug){
    std::cout << "easting " << easting << " northing " << northing << " e_coord " << e_coord << " n_coord " << n_coord << " x_min " << x_min << " y_min " << y_min << " cell size " << cell_size << std::endl;
  }
}



//_______________________________________________________________________________

void RampdemReader::LonLatToEastingNorthing(Double_t lon,Double_t lat,Double_t &easting,Double_t &northing){

  Double_t lon_rad = lon * TMath::DegToRad(); //convert to radians
  Double_t lat_rad = -lat * TMath::DegToRad();

  R_factor = scale_factor*c_0 * pow(( (1 + eccentricity*sin(lat_rad)) / (1 - eccentricity*sin(lat_rad)) ),eccentricity/2) * tan((TMath::Pi()/4) - lat_rad/2);

  easting = R_factor * sin(lon_rad);///(x_max-x_min);
  northing = R_factor * cos(lon_rad);///(y_max-y_min);

}



//_______________________________________________________________________________
void RampdemReader::ENtoLonLat(Int_t e_coord, Int_t n_coord, Double_t& lon, Double_t& lat) {
  //Takes as input the indicies from a BEDMAP data set, and turns them into latitude and longitude coordinates.  Information on which data set (surface data, ice depth, water depth) is necessary, in the form of coordinates of a corner of the map.  Code by Stephen Hoover.

  Double_t isometric_lat=0;
  Double_t easting = x_min+(cell_size*(e_coord+0.5)); //Add offset of 0.5 to get coordinates of middle of cell instead of edges.
  Double_t northing = -1*(y_min+(cell_size*(n_coord+0.5)));

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
    R_factor = TMath::Abs(easting/sin(lon));
  else if (easting == 0 && northing != 0)
    R_factor = TMath::Abs(northing);
  else {
    lat = 0; //at the pole, set lat=0 degrees
    lon = lon*TMath::RadToDeg();
    return;
  } //else

  isometric_lat = (TMath::Pi()/2) - 2*atan(R_factor/(scale_factor*c_0));

  lat = isometric_lat + a_bar*sin(2*isometric_lat) + b_bar*sin(4*isometric_lat) + c_bar*sin(6*isometric_lat) + d_bar*sin(8*isometric_lat);

  lon = lon * TMath::RadToDeg();  //convert to degrees
  lat =  -lat*TMath::RadToDeg(); //convert to degrees, with -90 degrees at the south pole
  return;
} //method ENtoLonLat



//_______________________________________________________________________________
void RampdemReader::EastingNorthingToLonLat(Double_t easting,Double_t northing,Double_t &lon,Double_t &lat){

  Int_t e_coord;
  Int_t n_coord;

  EastingNorthingToEN(easting,northing,e_coord,n_coord);
  ENtoLonLat(e_coord,n_coord,lon,lat);

  return;

}







//_______________________________________________________________________________

TProfile2D *RampdemReader::rampMap(int coarseness_factor, int set_log_scale, UInt_t &xBins, UInt_t &yBins){

  UInt_t num_columns = surface_elevation.size()/coarseness_factor;
  UInt_t num_rows = surface_elevation[0].size()/coarseness_factor;
  xBins = num_columns;
  yBins = num_rows;

  std::cout << "xBins " << xBins << " yBins " << yBins << std::endl;
  std::cout << "x_min " << x_min << "\t"
	    << "x_max " << x_max << "\t"
	    << "y_min " << y_min << "\t"
	    << "y_max " << y_max << "\t"
	    << std::endl;

  TProfile2D * theHist = new TProfile2D("h_antarctica_surface_elevation_hist", "",
					num_columns, x_min, x_max,
					num_rows, y_min, y_max+cell_size);

  for (unsigned int y=0; y < surface_elevation[0].size(); ++y){
    for (unsigned int x=0; x < surface_elevation.size(); ++x){
      if (set_log_scale == 1){
	theHist->Fill( x_min + double(x)*cell_size, -(y_min + double(y)*cell_size), log10(surface_elevation[x][y]+1000) ); //The "+1000" makes sure everything is positive.
      }
      else{
	theHist->Fill( x_min + double(x)*cell_size, -(y_min + double(y)*cell_size), surface_elevation[x][y] );
      }
    } //end for (loop over surface elevation vector)
  }

  theHist->SetStats(0);

  return theHist;

}




TProfile2D *RampdemReader::rampMapPartial(int coarseness_factor,double centralLon,double centralLat,double rangeMetres,Int_t &xBins,Int_t &yBins,Double_t &xMin,Double_t &xMax,Double_t &yMin,Double_t &yMax){

  Bool_t debug=false;

//   TColor mapColors;
//   Int_t mapColorInts[20];
//   char colorName[FILENAME_MAX];
//   for(int i=0;i<20;i++){
//     sprintf(colorName,"mapColor_%d",i);
//     if(i==0){
//       //mapColors.SetRGB(0,0,255);
//       mapColorInts[i] = mapColors.GetColor(51,51,255);
//     }
//     else{
//       //mapColors.SetRGB(204-60+i*3,255,255);
//       mapColorInts[i] = mapColors.GetColor(255-133+i*7,255,255);
//     }
//     //mapColorInts[i] = 10+i;
//   }

//   gStyle->SetPalette(20,mapColorInts);

  Int_t central_e_coord,central_n_coord;
  Int_t max_e_coord,min_e_coord;
  Int_t max_n_coord,min_n_coord;
  Double_t central_easting,central_northing;
  LonLattoEN(centralLon,centralLat,central_e_coord,central_n_coord);
  LonLatToEastingNorthing(centralLon,centralLat,central_easting,central_northing);
  EastingNorthingToEN(central_easting+rangeMetres,central_northing+rangeMetres,max_e_coord,max_n_coord);
  EastingNorthingToEN(central_easting-rangeMetres,central_northing-rangeMetres,min_e_coord,min_n_coord);


  if(debug){
    std::cout << "e_coord: min " << min_e_coord << " central " << central_e_coord << " max " << max_e_coord << std::endl;
    std::cout << "n_coord: min " << min_n_coord << " central " << central_n_coord << " max " << max_n_coord << std::endl;
  }

  if(max_e_coord<min_e_coord){
    Int_t swap_e = max_e_coord;
    max_e_coord = min_e_coord;
    min_e_coord = swap_e;
  }
  if(max_n_coord<min_n_coord){
    Int_t swap_n = max_n_coord;
    max_n_coord = min_n_coord;
    min_n_coord = swap_n;
  }
  xBins = max_e_coord-min_e_coord;
  yBins = max_n_coord-min_n_coord;
  xMin = central_easting-rangeMetres;
  xMax = central_easting+rangeMetres;
  yMin = central_northing-rangeMetres;
  yMax = central_northing+rangeMetres;

  char histName[FILENAME_MAX];
  sprintf(histName,"antarctica_surface_elevation_partial");
  TProfile2D *theHist = new TProfile2D(histName,"",(max_e_coord-min_e_coord)/coarseness_factor,central_easting-rangeMetres,central_easting+rangeMetres,(max_n_coord-min_n_coord)/coarseness_factor,central_northing-rangeMetres,central_northing+rangeMetres);

  if(debug){
    std::cout << "e_coord: min " << min_e_coord << " central " << central_e_coord << " max " << max_e_coord << std::endl;
    std::cout << "n_coord: min " << min_n_coord << " central " << central_n_coord << " max " << max_n_coord << std::endl;
  }

  for(Int_t row_index=min_n_coord;row_index<=max_n_coord;++row_index){
    if(row_index>=(int)(surface_elevation[0].size())) continue;
    for(Int_t column_index=min_e_coord;column_index<=max_e_coord;++column_index){

//       if(debug){
// 	std::cout << "row " << row_index << " up to " << max_n_coord << std::endl;
// 	std::cout << "col " << column_index << " up to " << max_e_coord << std::endl;
//       }

      if(column_index>=(int)(surface_elevation.size())) continue;

      theHist->Fill( x_min + double(column_index)*cell_size, -(y_min + double(row_index)*cell_size), surface_elevation[column_index][row_index]);
    }
  }

  theHist->SetStats(0);
  return theHist;
}




TGaxis *RampdemReader::distanceScale(Double_t xMin,Double_t xMax,Double_t yMin,Double_t yMax){
  //here xmin, xmax etc are the positions in easting/northing of the distance scale
  if(xMax<=xMin && yMax<=yMin) {
    std::cerr << "size ordering wrong: xMin " << xMin << " xMax " << xMax << " yMin " << yMin << " yMax " << yMax << std::endl;
    return NULL;
  }

  TGaxis *theAxis = new TGaxis(xMin,yMin,xMax,yMax,0,sqrt((xMax-xMin)/1e3*(xMax-xMin)/1e3-(yMax-yMin)/1e3*(yMax-yMin)/1e3),2,"");
  theAxis->SetTitle("km");

  return theAxis;
}




void RampdemReader::getMapCoordinates(double &xMin,double &yMin,double &xMax,double &yMax){

  xMin = x_min;
  yMin = y_min;
  xMax = x_max;
  yMax = y_max+cell_size;

}










////////////////////////////////////////////////////////////////////////////////////
//
//                                   BEDMAP 2 STUFF
//
////////////////////////////////////////////////////////////////////////////////////





static const char* dataSetToString(RampdemReader::dataSet dataSet){

  switch(dataSet){
  case RampdemReader::rampdem:
    return "rampdem";
  case RampdemReader::bed:
    return "bed";
  case RampdemReader::coverage:
    return "coverage";
  case RampdemReader::grounded_bed_uncertainty:
    return "grounded_bed_uncertainty";
  case RampdemReader::icemask_grounded_and_shelves:
    return "icemask_grounded_and_shelves";
  case RampdemReader::lakemask_vostok:
    return "lakemask_vostok";
  case RampdemReader::rockmask:
    return "rockmask";
  case RampdemReader::surface:
    return "surface";
  case RampdemReader::thickness:
    return "thickness";
  case RampdemReader::bedmap2_thickness_uncertainty_5km:
    return "bedmap2_thickness_uncertainty_5km";
  default:
    std::cerr << "Error in " << __FILE__ << ", unknown RampdemReader::dataSet requested" << std::endl;
    return NULL;

  }
}




static const char* dataSetToAxisTitle(RampdemReader::dataSet dataSet){
  // from the bedmap2 readme

  switch(dataSet){
  case RampdemReader::rampdem:
    return "Surface height (m)";
  case RampdemReader::bed:
    return "Bed height (m)";
  case RampdemReader::coverage:
    return "Ice coverage data";
  case RampdemReader::grounded_bed_uncertainty:
    return "Bed Uncertainty Grid (m)";
  case RampdemReader::icemask_grounded_and_shelves:
    return "Grounding line and floating ice shelves";
  case RampdemReader::lakemask_vostok:
    return "";
  case RampdemReader::rockmask:
    return "Rock Outcrops (m)";
  case RampdemReader::surface:
    return "Surface height (m)";
  case RampdemReader::thickness:
    return "Ice thickness (m)";
  case RampdemReader::bedmap2_thickness_uncertainty_5km:
    return "Ice thickness uncertainty";
  default:
    std::cerr << "Error in " << __FILE__ << ", unknown RampdemReader::dataSet requested" << std::endl;
    return NULL;

  }
}




static const VecVec& getDataIfNeeded(RampdemReader::dataSet dataSet){

  // If we haven't initialized the map with empty vectors, do it here
  if(bedMap2Data.size()==0){
    bedMap2Data[RampdemReader::rampdem] = VecVec();
    bedMap2Data[RampdemReader::bed] = VecVec();
    bedMap2Data[RampdemReader::coverage] = std::vector<std::vector<short> >();
    bedMap2Data[RampdemReader::grounded_bed_uncertainty] = std::vector<std::vector<short> >();
    bedMap2Data[RampdemReader::icemask_grounded_and_shelves] = std::vector<std::vector<short> >();
    bedMap2Data[RampdemReader::rockmask] = std::vector<std::vector<short> >();
    bedMap2Data[RampdemReader::surface] = std::vector<std::vector<short> >();
    bedMap2Data[RampdemReader::thickness] = std::vector<std::vector<short> >();
  }

  DataMap::iterator it = bedMap2Data.find(dataSet);
  if(it==bedMap2Data.end()){
    std::cerr << "Error in " << __FILE__ << ", unable to find requested data set!" << std::endl;
    // std::cerr << it->first << std::endl;
  }

  std::vector<std::vector<short >>& data = it->second;

  if(data.size() == 0){
    const char* anitaEnv = "ANITA_UTIL_INSTALL_DIR";
    const char* anitaUtilInstallDir = getenv(anitaEnv);
    if(anitaUtilInstallDir==NULL){
      std::cerr << "Error in " << __FILE__ << ", could not find environment variable " << anitaEnv << std::endl;
    }


    // Start with the anita install directory...
    std::string fileName(anitaUtilInstallDir);
    // ... append the calib subdir
    fileName.append("/share/anitaCalib/bedmap2_bin/bedmap2_");

    // append the appropriate filename
    const char* dataName = dataSetToString(dataSet);
    fileName.append(dataName);

    // Header file suffix
    std::string headName = fileName;
    headName.append(".hdr");


    // Open header file
    std::ifstream header(headName);
    if(!header.is_open()){
      std::cerr << "Error! Unable to open file " << headName << std::endl;
    }
    else{

      // Parse and store variables
      do{
	std::string key, value;
	header >> key >> value;

	// Think these are all we care about...
	if(key=="ncols"){
	  numXs[dataSet] = atoi(value.c_str());
	}
	else if(key=="nrows"){
	  numYs[dataSet] = atoi(value.c_str());
	}
	else if(key=="NODATA_value"){
	  noDatas[dataSet] = atoi(value.c_str());
	}
	else if(key=="xllcorner"){
	  minXs[dataSet]= atoi(value.c_str());
	}
	else if(key=="yllcorner"){
	  minYs[dataSet] = atoi(value.c_str());
	}
	else if(key=="cellsize"){
	  cellSizes[dataSet] = atoi(value.c_str());
	}
	std::cout << key << "\t" << value << "\t" << header.eof() << std::endl;
      } while(!header.eof());


      // Calculate other edge points.
      maxXs[dataSet] = minXs[dataSet] + numXs[dataSet]*cellSizes[dataSet];
      maxYs[dataSet] = minYs[dataSet] + numYs[dataSet]*cellSizes[dataSet];

      // Now get data file...
      fileName.append(".flt");
      FILE* fBedMap2 = fopen(fileName.c_str(), "r");
      if(fBedMap2==NULL){
	std::cerr << "Error in " << __FILE__ << ", could not open file " << fileName << std::endl;
      }
      else{
	// int nCols = bedMap2Headers[HeaderKey(dataSet, "nrows")];
	// int nRows = bedMap2Headers[HeaderKey(dataSet, "ncols")];

	const int numX = numXs[dataSet];
	const int numY = numYs[dataSet];

	std::vector<float> tempData(numX, 0);

	for(int y=0; y < numY; y++){
	  fread(&tempData[0], sizeof(float), numX, fBedMap2);
	  data.push_back(std::vector<short>(numX, 0));

	  for(int x = 0; x < numX; x++){
	    // std::cout << y << "\t" << x << "\t" << data.size() << "\t" << data.at(y).size() << "\t" << std::endl;
	    data.at(y).at(x) = short(tempData.at(x));
	  }
	}

	fclose(fBedMap2);
      }
    }
  }
  return data;
}




TH2D* RampdemReader::getHist(RampdemReader::dataSet dataSet, int coarseness){

  const VecVec& data = getDataIfNeeded(dataSet);

  Double_t cellSize = cellSizes[dataSet];

  Int_t xBins = numXs[dataSet]/coarseness;
  Int_t yBins = numYs[dataSet]/coarseness;

  Int_t xMin = minXs[dataSet];
  Int_t yMin = minYs[dataSet];


  Double_t yMax = yMin + numYs[dataSet]*cellSize;
  Double_t xMax = xMin + numXs[dataSet]*cellSize;

  Int_t noData = noDatas[dataSet];

  std::cout << "xMin " << xMin << "\t"
	    << "xMax " << xMax << "\t"
	    << "yMin " << yMin << "\t"
	    << "yMax " << yMax << "\t"
	    << std::endl;

  // VecVec& data = bedMap2Data.find(dataSet);


  // UInt_t num_columns = surface_elevation.size()/coarseness_factor;
  // UInt_t num_rows = surface_elevation[0].size()/coarseness_factor;
  // xBins = num_columns;
  // yBins = num_rows;

  // std::cout << "xBins " << xBins << " yBins " << yBins << std::endl;

  TString hName = TString::Format("h_%s_%d", dataSetToString(dataSet), coarseness);
  TProfile2D * theHist = new TProfile2D(hName, "", xBins, xMin, xMax, yBins, yMin, yMax);

  for(Int_t yBin=0; yBin < yBins; yBin++){
    for(Int_t xBin=0; xBin < xBins; xBin++){
      Double_t meanVal = 0;
      Int_t numVals = 0;

      for(int j=yBin*coarseness; j < (yBin+1)*coarseness; j++){
	for(int i=xBin*coarseness; i < (xBin+1)*coarseness; i++){


	  if(data[j][i] != noData){

	    meanVal += data[j][i];
	    numVals++;
	    // std::cout << j << "\t" << i << "\t" << data[j][i] << "\t" << numVals << std::endl;
	  }
	}
      }
      if(numVals > 0){
	theHist->Fill(xMin + double(xBin)*cellSize*coarseness, -(yMin + double(yBin)*cellSize*coarseness),
			meanVal/numVals);
	// theHist->SetBinContent(xBin + 1, yBins - yBin,
	// 		       // meanVal/numVals);
	// 		       10);

      }
    }
  }

  theHist->SetStats(0);

  return theHist;

}
