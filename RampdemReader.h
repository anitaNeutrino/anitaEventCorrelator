//////////////////////////////////////////////////////////////////////////////
/////  RampdemReader.h       Rampdem Data Reader                         /////
/////                                                                    /////
/////  Description:                                                      /////
/////     Some functions to get surface elevation required for event     /////
/////     Almost all the code stolen from Stephen ...                    /////
/////                                                                    /////
/////  Author: Matt Mottram (mottram@hep.ucl.ac.uk)                      /////
//////////////////////////////////////////////////////////////////////////////

#ifndef RAMPDEMREADER_H
#define RAMPDEMREADER_H

//Includes
#include <iostream>
#include "TObject.h"
#include "TMath.h"
#include "TProfile2D.h"
#include <vector>


class RampdemReader// : public TObject
{

 public:

  RampdemReader();
  ~RampdemReader();

  static RampdemReader*  Instance(); ///<Instance generator


  Double_t Geoid(Double_t latitude);

  //BEDMAP data input methods
  int readRAMPDEM();


  //BEDMAP utility methods
  Double_t Area(Double_t latitude);

  void ENtoLonLat(Int_t e_coord, 
		  Int_t n_coord,
		  Double_t& lon, 
		  Double_t& lat);
  void LonLattoEN(Double_t lon, 
		  Double_t lat,
		  int& e_coord, 
		  int& n_coord);
  void EastingNorthingToEN(Double_t easting,
			   Double_t northing,
			   Int_t &e_coord,
			   Int_t &n_coord);
  void LonLatToEastingNorthing(Double_t lon,
			       Double_t lat,
			       Double_t &easting,
			       Double_t &northing);
  
  //Data Output methods
  
  Double_t Surface(Double_t longitude, Double_t latitude);
  Double_t SurfaceAboveGeoid(Double_t longitude, Double_t latitude);
  Double_t SurfaceAboveGeoidRampDem(Double_t longitude, Double_t latitude);

  TProfile2D *rampMap(int coarseness_factor, int set_log_scale,UInt_t &xBins,UInt_t &yBins);
  TProfile2D *rampMapPartial(int coarseness_factor,double centralLon,double centralLat,double rangeMetres,Int_t &xBins,Int_t &yBins);

  //Generic method to flip Endianness.
  //WARNING: Flips byte order of anything put in - do not use on things like stuctures or classes!
  template <class thing>
    inline void flipEndian(thing &in) 
    {
      int size = sizeof(thing);

      thing out;

      char* p_in = (char *) &in;
      char* p_out = (char *) &out;

      for(int i=0;i<size;i++) 
	p_out[i] = p_in[size-1-i];

      in = out;

      return;
    } //template <class thing> inline void AnalysisTools::flipEndian(thing &in) 



 protected:
  static RampdemReader *fgInstance;  
   // protect against multiple instances

 private:

  /** RAMP DEM data.  Note: x increases to the right, y increases downward.  **/
  std::vector< std::vector<short> > surface_elevation;
  double cell_size;
  double x_min;
  double x_max;
  double y_min;
  double y_max;
  int nRows_surface;
  int  nCols_surface;
  int  nBytes_surface;




  static double scale_factor;  //scale factor at pole corresponding to 71 deg S latitude of true scale (used in both BEDMAP and the RAMP DEM)
  static double ellipsoid_inv_f; //of Earth
  static double ellipsoid_b;
  static double eccentricity;
  static double a_bar;
  static double b_bar;
  static double c_bar;
  static double d_bar;
  static double c_0;
  static double R_factor;
  static double nu_factor;

};

#endif //RAMPDEMREADER_H

