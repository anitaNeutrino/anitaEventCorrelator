#ifndef GEOMAGNETIC_H
#define GEOMAGNETIC_H

#include <vector>
#include "TArrow.h"
#include "TVector3.h"

/** 
 * @namespace Functions to calculate the Earth's geo-magnetic field for ANITA
 *
 * Currently uses the IGRF model
 * 
 */

class TCanvas;
namespace GeoMagnetic{


/** 
 * @class to hold the results of geomagnetic field calculations
 *
 * Since the field components are in a spherical coordinate system, the unit vectors
 * move as a function of position. So here I keep the field result with the position.
 * 
 */
class FieldPoint : public TArrow {
 public:
  FieldPoint (UInt_t unixTime, double lon, double lat, double alt);
  virtual ~FieldPoint(){}
  virtual void Draw(Option_t* opt = "");

  double posX() const {return fPosition.X();}
  double posY() const {return fPosition.Y();}
  double posZ() const {return fPosition.Z();}
  double posR() const {return fPosition.Mag();}
  double posTheta() const {return fPosition.Theta();}
  double posPhi() const {return fPosition.Phi();}
  
  double componentX() const {return fField.X();}
  double componentY() const {return fField.Y();}
  double componentZ() const {return fField.Z();}

  // TODO, maybe for checking
  // double componentRHat() const;
  // double componentThetaHat() const;
  // double componentPhiHat() const;   
  
 private:
  double fDrawScaleFactor; ///< Conversion factor (metres/nT) to draw on an AntarcticaBackground
  TVector3 fField; ///< Cartesian components of the geomagnetic field
  TVector3 fPosition; ///< Location of the magnetic field
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

  TCanvas* plotFieldAtAltitude(UInt_t unixTime, double altitude);

}
#endif
