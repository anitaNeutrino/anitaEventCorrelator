#ifndef GEOMAGNETIC_H
#define GEOMAGNETIC_H

#include <vector>
#include "TVector3.h"


/** 
 * @namespace Functions to calculate the Earth's geo-magnetic field for ANITA
 *
 * Currently uses the IGRF model
 * 
 */
namespace GeoMagnetic{


/** 
 * @class to hold the results of geomagnetic field calculations
 *
 * Since the field components are in a spherical coordinate system, the unit vectors
 * move as a function of position. So here I keep the field result with the position.
 * 
 */
class Field {
 public:
  /** 
   * @brief Default constructor
   */
  Field () : X(0), Y(0), Z(0), r(0), theta(0), phi(0) {}
  /** 
   * Default destructor
   */
  virtual ~Field(){}
  double X; ///< northwards pointing magnetic field in nT
  double Y; ///< East/west pointing magnetic field in nT
  double Z; ///< Downwards pointing magnetic field in nT
  double r; ///< Radial position of vector origin in Earth centred spherical polar coordinates
  double theta; ///< Elevation angle (0 = along +ve Z axis to north pole) of vector origin
  double phi; ///< Azimuth angle (0 = Greenwich maridian), increases +ve is west.
};


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

  Field getFieldAtLonLatAlt(UInt_t unixTime, double lon, double lat, double alt);
  Field getFieldAtSpherical(UInt_t unixTime, double r, double theta, double phi);

}
#endif
