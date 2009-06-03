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


class BedmapReader// : public TObject
{

 public:

  BedmapReader();
  ~BedmapReader();

  static BedmapReader*  Instance(); ///<Instance generator

  //BEDMAP data
  Double_t surface_elevation[1096][911]; //elevation of surface above geoid


  Double_t Geoid(Double_t latitude);

  //BEDMAP data input methods
  void ReadSurfaceElevation();
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
  
 protected:
  static BedmapReader *fgInstance;  
   // protect against multiple instances

 private:


};

#endif //BEDMAPREADER_H
