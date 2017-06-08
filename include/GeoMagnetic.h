#ifndef GEOMAGNETIC_H
#define GEOMAGNETIC_H

#include <vector>
#include "TVector3.h"

namespace GeoMagnetic{




  double g(UInt_t unixTime, int n, int m);
  double h(UInt_t unixTime, int n, int m);

  double getPotentialAtSpherical(UInt_t unixTime, double r, double theta, double phi);
  double getPotentialAtLonLatAlt(UInt_t unixTime, double lon, double lat, double alt);

  double X_atLonLatAlt(UInt_t unixTime, double lon,  double lat, double alt);
  double X_atSpherical(UInt_t unixTime, double r,  double theta, double phi);

  double Y_atLonLatAlt(UInt_t unixTime, double lon,  double lat, double alt);
  double Y_atSpherical(UInt_t unixTime, double r,  double theta, double phi);

  double Z_atLonLatAlt(UInt_t unixTime, double lon,  double lat, double alt);
  double Z_atSpherical(UInt_t unixTime, double r,  double theta, double phi);

}
#endif
