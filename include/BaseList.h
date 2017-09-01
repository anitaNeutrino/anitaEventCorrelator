#ifndef BASELIST_H
#define BASELIST_H

#include "TString.h"
#include <vector>
#include "RampdemReader.h" 

namespace BaseList{

  class base{
  public:
    base(const TString& theName, const TString& source, double lat, double lon, double alt=0){
      name = theName;
      dataSource = source;
      latitude = lat;
      longitude = lon;
      altitude = alt;
      RampdemReader::LonLatToEastingNorthing(lon,lat,Easting,Northing); 

    }
    base(const TString& theName, double lat, double lon, double alt=0){
      name = theName;
      latitude = lat;
      longitude = lon;
      altitude = alt;
      RampdemReader::LonLatToEastingNorthing(lon,lat,Easting,Northing); 
    }

    TString name;
    TString dataSource;
    double latitude;
    double longitude;
    double altitude;
    double Easting; 
    double Northing; 
  };

  const base& getBase(UInt_t i);
  void makeBaseList(); 
  void makeEmptyBaseList(); 
  
  size_t getNumBases();

};


#endif // BASELIST_H
