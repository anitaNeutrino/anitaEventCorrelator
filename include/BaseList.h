#ifndef BASELIST_H
#define BASELIST_H

#include "TString.h"
#include <vector>

namespace BaseList{

  class base{
  public:
    base(const TString& theName, const TString& source, double lat, double lon, double alt=0){
      name = theName;
      dataSource = source;
      latitude = lat;
      longitude = lon;
      altitude = alt;
    }
    base(const TString& theName, double lat, double lon, double alt=0){
      name = theName;
      latitude = lat;
      longitude = lon;
      altitude = alt;
    }

    TString name;
    TString dataSource;
    double latitude;
    double longitude;
    double altitude;
  };

  const base& getBase(UInt_t i);
  void makeBaseList();
  void makeEmptyBaseList();
  size_t getNumBases();

};


#endif // BASELIST_H
