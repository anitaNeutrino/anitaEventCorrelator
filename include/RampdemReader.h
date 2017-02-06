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
#include "TGaxis.h"
#include <vector>
#include <map>





class RampdemReader// : public TObject
{

 public:


  // Enumerates the sets of rampdem/bedmap2 data
  typedef enum{
    rampdem,
    bed,
    coverage,
    grounded_bed_uncertainty,
    icemask_grounded_and_shelves,
    lakemask_vostok,
    rockmask,
    surface,
    thickness,
    bedmap2_thickness_uncertainty_5km
  } dataSet;


  RampdemReader();
  ~RampdemReader();

  static RampdemReader*  Instance(); ///<Instance generator


  static Double_t Geoid(Double_t latitude);

  //BEDMAP data input methods
  static int readRAMPDEM();


  //BEDMAP utility methods
  static Double_t Area(Double_t latitude);

  static void ENtoLonLat(Int_t e_coord,
			 Int_t n_coord,
			 Double_t& lon,
			 Double_t& lat);
  static void LonLattoEN(Double_t lon,
			 Double_t lat,
			 int& e_coord,
			 int& n_coord);
  static void EastingNorthingToEN(Double_t easting,
				  Double_t northing,
				  Int_t &e_coord,
				  Int_t &n_coord);
  static void LonLatToEastingNorthing(Double_t lon,
				      Double_t lat,
				      Double_t &easting,
				      Double_t &northing);
  static void EastingNorthingToLonLat(Double_t easting,
				      Double_t northing,
				      Double_t &lon,
				      Double_t &lat);

  //Data Output methods

  static Double_t Surface(Double_t longitude, Double_t latitude);
  static Double_t SurfaceAboveGeoid(Double_t longitude, Double_t latitude);
  static Double_t SurfaceAboveGeoidRampDem(Double_t longitude, Double_t latitude);

  static TProfile2D *rampMap(int coarseness_factor, int set_log_scale,UInt_t &xBins,UInt_t &yBins);
  static TProfile2D *rampMapPartial(int coarseness_factor,double centralLon,double centralLat,double rangeMetres,Int_t &xBins,Int_t &yBins,Double_t &xMin,Double_t &xMax,Double_t &yMin,Double_t &yMax);
  static TGaxis *distanceScale(Double_t xMin,Double_t xMax,Double_t yMin,Double_t yMax);

  //Generic method to flip Endianness.
  //WARNING: Flips byte order of anything put in - do not use on things like stuctures or classes!
  template <class thing>
  inline static void flipEndian(thing &in)
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


  static void getMapCoordinates(double &xMin,double &yMin,double &xMax,double &yMax);
  static TH2D* getHist(RampdemReader::dataSet dataSet, int coarseness_factor);

 protected:
  static RampdemReader *fgInstance;
   // protect against multiple instances


};

#endif //RAMPDEMREADER_H
