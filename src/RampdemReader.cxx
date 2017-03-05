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



// Typedefs for parsing the surface data
typedef std::vector<std::vector<short> > VecVec;
typedef std::map<RampdemReader::dataSet, VecVec > DataMap;
static DataMap bedMap2Data;
typedef std::map<RampdemReader::dataSet, Double_t> HeaderMap;

static HeaderMap numXs;
static HeaderMap numYs;
static HeaderMap noDatas;
static HeaderMap minXs;
static HeaderMap minYs;
static HeaderMap maxXs;
static HeaderMap maxYs;
static HeaderMap cellSizes;


// static functions to read in the data / generic fill histogram function.
static const VecVec& getDataIfNeeded(RampdemReader::dataSet dataSet);
static TProfile2D* fillThisHist(TProfile2D* theHist, RampdemReader::dataSet dataSet);




//Variables for conversion between polar stereographic coordinates and lat/lon.
// Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf

// scale factor at pole corresponding to 71 deg S latitude of true scale (used in both BEDMAP and RAMP DEM)
static double scale_factor=0.97276901289;

static double ellipsoid_inv_f = 298.257223563; //of Earth

// static double ellipsoid_b = R_EARTH*(1-(1/ellipsoid_inv_f)); // Unused.

static double eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));

static double a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;

static double b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;

static double c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;

static double d_bar = 4279*pow(eccentricity,8)/161280;

static double c_0 = (2*R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);

// Varies with latitude, defined here for 71 deg S...
static double R_factor = scale_factor*c_0 * pow(( (1 + eccentricity*sin(71*TMath::RadToDeg())) / (1 - eccentricity*sin(71*TMath::RadToDeg())) ),eccentricity/2) * tan((TMath::Pi()/4) - (71*TMath::RadToDeg())/2);

static double nu_factor = R_factor / cos(71*TMath::RadToDeg());




RampdemReader*  RampdemReader::fgInstance = 0; //!< Pointer to instance.


















/**
 * Default constructor.
 * This class has been converted so that it works entirely statically
 * You don't need to call this or create an instance...
 * Preserved for backward compatibility.
 *
 */
RampdemReader::RampdemReader(){
  fgInstance=this;
}


/**
 * Default destructor, not used.
 */

RampdemReader::~RampdemReader(){}



/**
 * Instance generated. Deprecated.
 *
 *
 * @return pointer to RampdemReader singleton
 */
RampdemReader*  RampdemReader::Instance(){
  //static function
  return (fgInstance) ? (RampdemReader*) fgInstance : new RampdemReader();
}












/**
 * Returns the height of the Antarctic surface above the centre of the Earth.
 *
 * @param lon is the longitude (degrees)
 * @param lat is the latitude (degrees)
 *
 * @return height of the Antarctic surface above the geoid.
 */
Double_t RampdemReader::Surface(Double_t lon,Double_t lat) {
  return (SurfaceAboveGeoid(lon,lat) + Geoid(lat));
}







/**
 * Returns the elevation above the geoid of the surface of the top of the ice (or bare ground if no ice)
 * in meters, at a location specified by a latitude and longitude (in degrees).
 * @param lon is the longiutde (degrees)
 * @param lat is the latitude (degrees)
 *
 * @return elevation above geoid in metres
 */
Double_t RampdemReader::SurfaceAboveGeoid(Double_t lon, Double_t lat, RampdemReader::dataSet dataSet) {

  getDataIfNeeded(dataSet);
  VecVec& surface_elevation = bedMap2Data[dataSet];

  Int_t nCols_surface = numXs[dataSet];
  Int_t nRows_surface = numYs[dataSet];

  Double_t surface=0;

  Int_t e_coord_surface=0;
  Int_t n_coord_surface=0;
  LonLattoEN(lon,lat,e_coord_surface,n_coord_surface, dataSet);

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






/**
 * Returns the height of the Earth surface geoid in metres.
 *
 * @param latitude in degrees
 *
 * @return the height of the geoid in metres
 */
Double_t RampdemReader::Geoid(Double_t latitude) {
  return (GEOID_MIN*GEOID_MAX/sqrt(pow(GEOID_MIN,2)-(pow(GEOID_MIN,2)-pow(GEOID_MAX,2))*pow(cos(latitude*TMath::DegToRad()),2)));
}






/**
 * Function to read in the original RAMPDEM data.
 * This function is called by getDataIfNeeded, as the RAMPDEM elevation data files have a different format from the BEDMAP2 data
 *
 * @return greater than zero if reading/parsing files fails for some reason.
 */
int RampdemReader::readRAMPDEM(){

  // std::cerr << __PRETTY_FUNCTION__ <<  std::endl;

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

  double cell_size, x_min, x_max, y_min, y_max, min_value, max_value, mean, std_deviation;
  int nRows_surface, nCols_surface, nBytes_surface;
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
  cellSizes[RampdemReader::rampdem] = cell_size;
  noDatas[RampdemReader::rampdem] = -9999; // by hand

  // emptry VecVec is now initially put in by the getDataIfNeeded function
  // bedMap2Data[RampdemReader::rampdem] = VecVec();

  VecVec& surface_elevation = bedMap2Data[RampdemReader::rampdem];

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
      }
      flipEndian(temp_data);
      surface_elevation[column_index][row_index] = temp_data;
    }
  }
  if(dem_data.read((char*)&temp_data, sizeof(short))){
    std::cerr << "[RampdemReader::readRAMPDEM] Error! I could read a value (" << temp_data
	      << ") from the data file after I thought that it was empty!.\n";
    return 2;
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
Double_t RampdemReader::Area(Double_t latitude, RampdemReader::dataSet dataSet) {

  getDataIfNeeded(dataSet);
  Double_t cell_size = cellSizes[dataSet];

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
 * @param dataSet is which data set to get the coordinate of
 */
// void RampdemReader::LonLattoEN(Double_t lon, Double_t lat, int& e_coord, int& n_coord) {
void RampdemReader::LonLattoEN(Double_t lon, Double_t lat, int& e_coord, int& n_coord, RampdemReader::dataSet dataSet) {

  getDataIfNeeded(dataSet);

  Double_t easting=0;
  Double_t northing=0;

  LonLatToEastingNorthing(lon,lat,easting,northing);
  EastingNorthingToEN(easting,northing,e_coord,n_coord, dataSet);

}






/**
 * Converts Easting/northing to RAMPDEM data indices
 *
 * @param easting
 * @param northing
 * @param e_coord is the easting coordinate index
 * @param n_coord is the the northing coordinate index
 * @param dataSet is the data set to get the e/n coordinate of
 */
void RampdemReader::EastingNorthingToEN(Double_t easting,Double_t northing,Int_t &e_coord,Int_t &n_coord, RampdemReader::dataSet dataSet){

  getDataIfNeeded(dataSet);

  Int_t x_min = minXs[dataSet];
  Int_t y_min = minYs[dataSet];
  Int_t cell_size = cellSizes[dataSet];

  e_coord = (int)((easting - x_min) / cell_size);
  n_coord = (int)((-1*northing - y_min) / cell_size);

  // std::cout << easting << "\t" << northing << "\t" << e_coord << "\t" << n_coord << std::endl;
}






/**
 * Convert longitude and latitude to easting and northing using the geoid model
 *
 * @param lon is the longitude in degrees
 * @param lat is the latitude in degrees
 * @param easting in meters
 * @param northing in meters
 */
void RampdemReader::LonLatToEastingNorthing(Double_t lon,Double_t lat,Double_t &easting,Double_t &northing){

  Double_t lon_rad = lon * TMath::DegToRad(); //convert to radians
  Double_t lat_rad = -lat * TMath::DegToRad();

  R_factor = scale_factor*c_0 * pow(( (1 + eccentricity*sin(lat_rad)) / (1 - eccentricity*sin(lat_rad)) ),eccentricity/2) * tan((TMath::Pi()/4) - lat_rad/2);

  easting = R_factor * sin(lon_rad);///(x_max-x_min);
  northing = R_factor * cos(lon_rad);///(y_max-y_min);

}





/**
 * Takes as input the indicies from a BEDMAP data set, and turns them into latitude and longitude coordinates.
 * Original code by Stephen Hoover.
 *
 * @param e_coord is the easting coordinate
 * @param n_coord is the northing coordinate
 * @param lon is the longiude in degrees
 * @param lat is the latitude in degrees
 * @param dataSet picks the data set to get the coordinates at the specified lat/lon for.
 */
void RampdemReader::ENtoLonLat(Int_t e_coord, Int_t n_coord, Double_t& lon, Double_t& lat, RampdemReader::dataSet dataSet) {
  //

  getDataIfNeeded(dataSet);
  Double_t x_min = minXs[dataSet];
  Double_t y_min = minYs[dataSet];
  Double_t cell_size = cellSizes[dataSet];

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






/**
 * Convert from easting/northing to longitude and latitude
 *
 * @param easting in meters
 * @param northing in meters
 * @param lon is the longitude
 * @param lat is the latitude
 */
void RampdemReader::EastingNorthingToLonLat(Double_t easting,Double_t northing,Double_t &lon,Double_t &lat, RampdemReader::dataSet dataSet){

  Int_t e_coord;
  Int_t n_coord;

  EastingNorthingToEN(easting,northing,e_coord,n_coord, dataSet);
  ENtoLonLat(e_coord,n_coord,lon,lat, dataSet);

  return;

}







/**
 * DEPRECATED.
 * Left for backward compatibility, prefer getMap(RampdemReader::rampdem...)
 *
 * @param coarseness downsamples the easting/northing bins
 * @param set_log_scale converts heights to log units
 * @param xBins is the number of bins on the x-axis
 * @param yBins is the number of bins on the y-axis
 *
 * @return a histogram of the surface elevation from the RAMPDEM data set.
 */
TProfile2D *RampdemReader::rampMap(int coarseness, int set_log_scale, UInt_t &xBins, UInt_t &yBins){

  TProfile2D* theHist = getMap(RampdemReader::rampdem, coarseness);

  xBins = theHist->GetNbinsX();
  yBins = theHist->GetNbinsY();

  if(set_log_scale){
    for(UInt_t bx = 1; bx <= xBins; bx++){
      for(UInt_t by = 1; by <= yBins; by++){
	double val = theHist->GetBinContent(bx, by);
	val += 1000; // avoid negative numbers
	theHist->SetBinContent(bx, by, TMath::Log10(val));
      }
    }
  }

  theHist->SetStats(0);

  return theHist;

}




/**
 * DEPRECATED.
 * Left for backward compatibility, prefer getMapPartial(RampdemReader::rampdem...)
 *
 * @param coarseness downsamples the easting/northing bins
 * @param centralLon is the latitude (degrees) on which to centre the histogram
 * @param centralLat is the longitude (degrees) on which to centre the historam
 * @param rangeMetres is the extent of the histogram edges from the centre point (metres)
 * @param xBins is the number of bins on the x-axis
 * @param yBins is the number of bins on the y-axis
 * @param xMin is the lower limit of the x-axis
 * @param xMax is the upper limit of the x-axis
 * @param yMin is the lower limit of the y-axis
 * @param yMax is the upper limit of the y-axis
 *
 * @return a histogram of the surface elevation from the RAMPDEM data set in the specified region.
 */
TProfile2D *RampdemReader::rampMapPartial(int coarseness,
					  double centralLon, double centralLat, double rangeMetres,
					  Int_t &xBins, Int_t &yBins,
					  Double_t &xMin, Double_t &xMax,
					  Double_t &yMin,Double_t &yMax){


  TProfile2D* theHist = getMapPartial(RampdemReader::rampdem, coarseness, centralLon, centralLat, rangeMetres);
  xBins = theHist->GetNbinsX();
  yBins = theHist->GetNbinsY();
  xMin = theHist->GetXaxis()->GetBinLowEdge(1);
  xMax = theHist->GetXaxis()->GetBinLowEdge(xBins+1);
  yMin = theHist->GetYaxis()->GetBinLowEdge(1);
  yMax = theHist->GetYaxis()->GetBinLowEdge(yBins+1);

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



/**
 * Get the corners of the map (for a given data set)
 *
 * @param xMin is the minimum easting (metres)
 * @param yMin is the minimum northing (metres)
 * @param xMax is the maximum easting (metres)
 * @param yMax is the maximum northing (metres)
 * @param dataSet is the data set from which to select the coordinates
 */
void RampdemReader::getMapCoordinates(double &xMin, double &yMin,
				      double &xMax, double &yMax,
				      RampdemReader::dataSet dataSet){



  getDataIfNeeded(dataSet);
  xMin = minXs[dataSet];
  yMin = minYs[dataSet];
  xMax = maxXs[dataSet];
  yMax = maxYs[dataSet]+cellSizes[dataSet];

}




/**
 * Get number of bins in X and Y for any data set
 *
 * @param numX
 * @param numY
 * @param dataSet
 */
void RampdemReader::getNumXY(Int_t& numX, Int_t&numY,
			     RampdemReader::dataSet dataSet){

  getDataIfNeeded(dataSet);
  numX = numXs[dataSet];
  numY = numYs[dataSet];

}










////////////////////////////////////////////////////////////////////////////////////
//
//                                   BEDMAP 2 STUFF
//
////////////////////////////////////////////////////////////////////////////////////




/**
 * Convert dataSet enum to string for reading in files.
 *
 * @param dataSet is the data set
 *
 * @return a c string containing the dataSet enum name.
 */
static const char* dataSetToString(RampdemReader::dataSet dataSet){

  switch(dataSet){
  case RampdemReader::rampdem:
    return "rampdem";
  case RampdemReader::bed:
    return "bed";
  // case RampdemReader::coverage:
  //   return "coverage";
  // case RampdemReader::grounded_bed_uncertainty:
  //   return "grounded_bed_uncertainty";
  case RampdemReader::icemask_grounded_and_shelves:
    return "icemask_grounded_and_shelves";
  // case RampdemReader::lakemask_vostok:
  //   return "lakemask_vostok";
  // case RampdemReader::rockmask:
  //   return "rockmask";
  case RampdemReader::surface:
    return "surface";
  case RampdemReader::thickness:
    return "thickness";
  // case RampdemReader::bedmap2_thickness_uncertainty_5km:
  //   return "bedmap2_thickness_uncertainty_5km";
  default:
    std::cerr << "Error in " << __FILE__ << ", unknown RampdemReader::dataSet requested" << std::endl;
    return NULL;
  }
}



/**
 * Get z-axis title for any data set
 *
 * @param dataSet is the data set
 *
 * @return the z-axis title
 */
const char* RampdemReader::dataSetToAxisTitle(RampdemReader::dataSet dataSet){
  // from the bedmap2 readme

  switch(dataSet){
  case RampdemReader::rampdem:
    return "Surface height (m)";
  case RampdemReader::bed:
    return "Bed height (m)";
  // case RampdemReader::coverage:
  //   return "Ice coverage data";
  // case RampdemReader::grounded_bed_uncertainty:
  //   return "Bed Uncertainty Grid (m)";
  case RampdemReader::icemask_grounded_and_shelves:
    return "Grounding line and floating ice shelves";
  // case RampdemReader::lakemask_vostok:
  //   return "lakemask_vostok";
  // case RampdemReader::rockmask:
  //   return "Rock Outcrops (m)";
  case RampdemReader::surface:
    return "Surface height (m)";
  case RampdemReader::thickness:
    return "Ice thickness (m)";
  // case RampdemReader::bedmap2_thickness_uncertainty_5km:
  //   return "Ice thickness uncertainty";
  default:
    std::cerr << "Error in " << __FILE__ << ", unknown RampdemReader::dataSet requested" << std::endl;
    return NULL;

  }
}





/**
 * Handles reading in any data set.
 * The BEDMAP2 data are handled in this function, RAMPDEM data read in by the old ReadRAMPDEM function
 *
 * @param dataSet is the selected data set
 *
 * @return a reference to the raw bedmap/rampdem data
 */

static const VecVec& getDataIfNeeded(RampdemReader::dataSet dataSet){


  // If we haven't initialized the map with empty vectors, do it here
  if(bedMap2Data.size()==0){
    bedMap2Data[RampdemReader::rampdem] = VecVec();
    bedMap2Data[RampdemReader::bed] = VecVec();
    // bedMap2Data[RampdemReader::coverage] = VecVec();
    // bedMap2Data[RampdemReader::grounded_bed_uncertainty] = VecVec();
    bedMap2Data[RampdemReader::icemask_grounded_and_shelves] = VecVec();
    // bedMap2Data[RampdemReader::lakemask_vostok] = VecVec();
    // bedMap2Data[RampdemReader::rockmask] = VecVec();
    bedMap2Data[RampdemReader::surface] = VecVec();
    bedMap2Data[RampdemReader::thickness] = VecVec();
  }


  DataMap::iterator it = bedMap2Data.find(dataSet);
  if(it==bedMap2Data.end()){
    std::cerr << "Error in " << __FILE__ << ", unable to find requested data set!" << std::endl;
    // std::cerr << it->first << std::endl;
  }

  VecVec& data = it->second;

  // special case for old rampdem data...
  if(dataSet==RampdemReader::rampdem){
    if(data.size()==0){
      RampdemReader::readRAMPDEM();
    }
    return data;
  }


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
	// std::cout << key << "\t" << value << "\t" << header.eof() << std::endl;
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
	    data.at(y).at(x) = short(tempData.at(x));
	  }
	}

	fclose(fBedMap2);
      }
    }
  }
  return data;
}






/**
 * Function to loop over a any created histogram and fill it with a given data set
 * This won't work unless you created it with the correct limits. (See
 *
 * @param theHist is the TProfile2D to be filled
 * @param dataSet is the selected data set
 *
 * @return the same histogram (theHist)
 */
TProfile2D* RampdemReader::fillThisHist(TProfile2D* theHist, RampdemReader::dataSet dataSet){

  const VecVec& data = getDataIfNeeded(dataSet);
  Double_t cellSize = cellSizes[dataSet];

  Int_t xBins = numXs[dataSet];
  Int_t yBins = numYs[dataSet];

  Double_t xMin = minXs[dataSet];
  Double_t yMin = minYs[dataSet];

  // Double_t xMax = maxXs[dataSet];
  // Double_t yMax = maxYs[dataSet];
  Double_t noData = noDatas[dataSet];

  if(dataSet==RampdemReader::rampdem){
    for(UInt_t yBin=0; yBin < data.at(0).size(); yBin++){
      for(UInt_t xBin=0; xBin < data.size(); xBin++){

	theHist->Fill(xMin + double(xBin)*cellSize,
		      -(yMin + double(yBin)*cellSize),
		      data.at(xBin).at(yBin));
      }
    }
  }
  else{
    for(Int_t yBin=0; yBin < yBins; yBin++){
      for(Int_t xBin=0; xBin < xBins; xBin++){
	if(data[yBin][xBin] != noData){
	  theHist->Fill(xMin + double(xBin)*cellSize,
			-(yMin + double(yBin)*cellSize),
			data[yBin][xBin]);
	}
      }
    }
  }

  theHist->SetStats(0);

  return theHist;

}




/**
 * Creates a lovely new map of Antarctica from the requested data set with the requested coarseness.
 *
 * @param dataSet is the dataSet
 * @param coarseness is a rebinning factor for X and Y
 *
 * @return the new TProfile2D
 */
TProfile2D* RampdemReader::getMap(RampdemReader::dataSet dataSet, int coarseness){

  getDataIfNeeded(dataSet);

  Int_t xBins = numXs[dataSet];
  Int_t yBins = numYs[dataSet];

  Double_t xMin = minXs[dataSet];
  Double_t yMin = minYs[dataSet];

  Double_t xMax = maxXs[dataSet];
  Double_t yMax = maxYs[dataSet];
  Double_t cellSize = cellSizes[dataSet];

  TString hName = TString::Format("h_%s_%d", dataSetToString(dataSet), coarseness);
  TProfile2D * theHist = new TProfile2D(hName, "",
					xBins/coarseness, xMin, xMax,
					yBins/coarseness, yMin, yMax+cellSize);


  fillThisHist(theHist, dataSet);
  return theHist;

}








/**
 * Create a TProfile2D of the continent with the specified data set, centred on the given coordinates.
 *
 * @param dataSet is the selected data set
 * @param coarseness downsamples the easting/northing bins
 * @param centralLon is the longitude (degrees) on which to centre the histogram
 * @param centralLat is the latitude (degrees) on which to centre the histogram
 * @param rangeMetres is the distance from the centre to any edge of the histogram
 *
 * @return the created histogram
 */
TProfile2D* RampdemReader::getMapPartial(RampdemReader::dataSet dataSet, int coarseness,
					 double centralLon, double centralLat, double rangeMetres){

  Int_t central_e_coord,central_n_coord;
  LonLattoEN(centralLon, centralLat, central_e_coord, central_n_coord, dataSet);

  Double_t central_easting,central_northing;
  LonLatToEastingNorthing(centralLon, centralLat, central_easting, central_northing);

  Int_t max_e_coord, max_n_coord;
  EastingNorthingToEN(central_easting + rangeMetres, central_northing + rangeMetres,
		      max_e_coord, max_n_coord, dataSet);

  Int_t min_e_coord, min_n_coord;
  EastingNorthingToEN(central_easting - rangeMetres, central_northing - rangeMetres,
		      min_e_coord, min_n_coord, dataSet);

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

  Int_t thisXBins = (max_e_coord-min_e_coord)/coarseness;
  Int_t thisYBins = (max_n_coord-min_n_coord)/coarseness;

  // std::cout << thisXBins << "\t" << thisYBins << std::endl;
  // std::cout << max_e_coord << "\t" << min_e_coord << std::endl;
  // std::cout << max_n_coord << "\t" << min_n_coord << std::endl;
  double thisXMin = central_easting - rangeMetres;
  double thisXMax = central_easting + rangeMetres;
  double thisYMin = central_northing - rangeMetres;
  double thisYMax = central_northing + rangeMetres;

  TString histName = TString::Format("h_%s_partial_%d", dataSetToString(dataSet), coarseness);
  TProfile2D *theHist = new TProfile2D(histName,"",
				       thisXBins, thisXMin, thisXMax,
				       thisYBins, thisYMin, thisYMax);

  fillThisHist(theHist, dataSet);

  return theHist;

}



/**
 * Returns true if a specified latitude or longitude is on the continent...
 * ...according to the RampdemReader::surface model
 *
 *
 * @param lon is the longitude (degrees)
 * @param lat is the latitude (degrees)
 *
 * @return true if on the continent, false otherwise
 */
Bool_t RampdemReader::isOnContinent(Double_t lon, Double_t lat){
  const VecVec& data = getDataIfNeeded(RampdemReader::surface);
  int e_coord, n_coord;
  LonLattoEN(lon, lat, e_coord, n_coord, RampdemReader::surface);
  // std::cout << lon << "\t" << lat << e_coord << "\t" << n_coord << std::endl;

  Bool_t isOnContinent = false;
  if(n_coord >= 0 && n_coord < numYs[RampdemReader::surface] &&
     e_coord >= 0 && e_coord < numXs[RampdemReader::surface]){

    if(data.at(n_coord).at(e_coord)!=noDatas[RampdemReader::surface]){
      isOnContinent = true;
    }

  }
  return isOnContinent;
}
