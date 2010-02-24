//////////////////////////////////////////////////////////////////////////////
/////  BedmapReader.h        Bedmap Data Reader                          /////
/////                                                                    /////
/////  Description:                                                      /////
/////     Some functions to get surface elevation required for event     /////
/////     relocation that don't really need to have their own class.     /////
/////                                                                    /////
/////  Author: Matt Mottram (mottram@hep.ucl.ac.uk)                      /////
//////////////////////////////////////////////////////////////////////////////

#ifndef BEDMAPREADER_H
#define BEDMAPREADER_H

//Includes
#include <iostream>
#include "TObject.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TGaxis.h"
#include <vector>



class BedmapReader// : public TObject
{

 public:

  BedmapReader(bool icethicknessMode);
  ~BedmapReader();

  static BedmapReader*  Instance(bool icethicknessMode); ///<Instance generator

  //BEDMAP data
/*   Double_t surface_elevation[1096][911]; //elevation of surface above geoid */
  std::vector< std::vector<short> > surface_elevation;
  //for now we're going to store ice thickness here - allows pretty maps like rampdemreader ones to be made, but with ice thickness instead of elevation - nice as it includes the ice shelfs


  Double_t Geoid(Double_t latitude);

  //BEDMAP data input methods
  void ReadSurfaceElevation(bool icethicknessMode);
  void ReadSurfaceElevationRampDem();


  //BEDMAP utility methods
  Double_t Area(Double_t latitude);

  void ENtoLonLat(Int_t e_coord, 
		  Int_t n_coord,
		  Double_t xLowerLeft,
		  Double_t yLowerLeft,
		  
		  Double_t& lon, 
		  Double_t& lat);
  void SurfaceENtoLonLat(Int_t e,
			 Int_t n,
			 
			 Double_t& lon,
			 Double_t& lat);
  void LonLattoEN(Double_t lon, 
		  Double_t lat,
		  Double_t xLowerLeft,
		  Double_t yLowerLeft,
		  
		  int& e_coord, 
		  int& n_coord);
  void SurfaceLonLattoEN(Double_t lon,
			 Double_t lat,
			 
			 int& e_coord,
			 int& n_coord);
  
  //Data Output methods
  
  Double_t Surface(Double_t longitude, Double_t latitude);
  Double_t SurfaceAboveGeoid(Double_t longitude, Double_t latitude);
  Double_t SurfaceAboveGeoidRampDem(Double_t longitude, Double_t latitude);
  TProfile2D *bedmapMap(int coarseness_factor, int set_log_scale,UInt_t &xBins,UInt_t &yBins);
  TGaxis *distanceScale(Double_t xMin,Double_t xMax,Double_t yMin,Double_t yMax);
  
 protected:
  static BedmapReader *fgInstance;  
   // protect against multiple instances

 private:

  Int_t cellSize;
  Int_t nCols_surface;
  Int_t nRows_surface;
  Int_t xLowerLeft_surface;
  Int_t yLowerLeft_surface;
  Int_t xUpperRight_surface;
  Int_t yUpperRight_surface;



};

#endif //BEDMAPREADER_H
