#ifndef GEOMAGNETIC_H
#define GEOMAGNETIC_H

#include <vector>
#include "TArrow.h"
#include "TVector3.h"
#include "UsefulAdu5Pat.h"

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
  FieldPoint (UInt_t unixTime, const TVector3& position);
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

  const TVector3& field(){return fField;}
  const TVector3& position(){return fPosition;}  
  UInt_t getUnixTime(){return fUnixTime;}
  
  // TODO, maybe for checking
  // double componentRHat() const;
  // double componentThetaHat() const;
  // double componentPhiHat() const;   
  
 private:
  void calculateFieldAtPosition();
  double fDrawScaleFactor; ///< Conversion factor (metres/nT) to draw on an AntarcticaBackground
  TVector3 fField; ///< Cartesian components of the geomagnetic field
  TVector3 fPosition; ///< Location of the magnetic field
  UInt_t fUnixTime;  
};


  double getExpectedPolarisation(UsefulAdu5Pat& usefulPat, double phiWave, double thetaWave);
  double getExpectedPolarisationUpgoing(UsefulAdu5Pat& usefulPat, double phiWave, double thetaWave, double pathLength);

  TVector3 getUnitVectorAlongThetaWavePhiWave(UsefulAdu5Pat& usefulPat, double phiWave, double thetaWave);

  const double n_air = 1;
  const double n_ice = 1.31; // might need to check this
  TVector3 specularReflection(const TVector3& reflectionPointToSource, const TVector3& surfaceNormal);
  TVector3 fresnelReflection(const TVector3& sourceToReflection, const TVector3& surfaceNormal, TVector3& electricFieldVec, double n1=n_air, double n2=n_ice);
  TCanvas* plotFresnelReflection();

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
  TCanvas* plotAtmosphere();

  double getAtmosphericDensity(double altitude);
  TVector3 getXMaxPosition(const TVector3& initialPosition, const TVector3& cosmicRayDirection);
  TVector3 getInitialPosition(const TVector3& destination, const TVector3& destinationToSource);
  void setDebug(bool db);


}
#endif
