//////////////////////////////////////////////////////////////////////////////
/////  UsefulAdu5Pat.h        Useful ANITA event class                      /////
/////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for utilising the ADU5 data                    /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#ifndef USEFULADU5PAT_H
#define USEFULADU5PAT_H

//Includes
#include <TObject.h>
#include <TGraph.h>
#include <TVector3.h>
#include "Adu5Pat.h"
#include "AnitaConventions.h"
#include "RampdemReader.h"
#include "AnitaGeomTool.h"

/**
 * @namespace AnitaStaticAdu5Offsets
 * @brief Static value to which the pitch/roll variables are set in UsefulAdu5Pat, currently for ANITA-3
 */
namespace AnitaStaticAdu5Offsets{
  const Double_t heading = 0;
  const Double_t pitch = 0;
  const Double_t roll = 0;
}

/**
 * @class UsefulAdu5Pat
 * @brief As UsefulAnitaEvent is to RawAnitaEvent, UsefulAdu5Pat is to Adu5Pat. Well not quite as useful but you get the picture.
 * 
 * This is the position and attitude class that inherits from Adu5Pat and has some useful methods for
 * determining the expected plane wave angles and inter-antenna crossing times for a given source.
 */

class UsefulAdu5Pat: public Adu5Pat
{

 public:
  UsefulAdu5Pat();			/// Default constructor
  UsefulAdu5Pat(const Adu5Pat *patPtr); /// Assignment constructor
  ~UsefulAdu5Pat();			/// Destructor

  /**
   * 
   * For a given azimuthal and elevation angle of a plane wave (in payload coordinates) calculates the point on the Earth's surface that the source would come from.
   * 
   * @param phiWave Azimuthal angle of plane wave (in payload centric coordinates)
   * @param thetaWave Elevation angle of plane wave (in payload centric coordinates)
   * @param sourceLon Reference to a Double_t in which the longitude of the source will be stored.
   * @param sourceLat Reference to a Double_t in which the latitude of the source will be stored.
   * @param desiredAlt Call specified altitude of the payload (if not specified, set to 0.0 so that this function would behave the same as getSourceLonAndLatAltZero() function)
   * @return 1 on successful source localisation, 0 if the solution does not point to the ground.
   */
  int getSourceLonAndLatAtDesiredAlt(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat, Double_t desiredAlt);

  
  int getSourceLonAndLatAltZero(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat);
  
  /**
   * 
   * For a given azimuthal and elevation angle of a plane wave (in payload coordinates) calculates the point on the Earth's surface that the source would come from
   * Accounts for the ice surface.
   * 
   * @param phiWave Azimuthal angle of plane wave (in payload centric coordinates)
   * @param thetaWave Elevation angle of plane wave (in payload centric coordinates)
   * @param sourceLon Reference to a Double_t in which the longitude of the source will be stored.
   * @param sourceLat Reference to a Double_t in which the latitude of the source will be stored.
   * @param desiredAlt Call specified altitude of the payload (if not specified, set to 0.0 so that this function would behave the same as getSourceLonAndLatAltZero() function)
   * @return 1 or 2 successful source localisation, -1, 0, or 4 for various failures, all of which set sourceLon/Lat/Alt to -9999 (see comments in source code for explanation)
   */
  int getSourceLonAndLatAtAlt(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat,Double_t &sourceAltitude);

  /**
   * Const version of getSourceLonAndLatAtAlt, all params and return values are the same.
   */
  int getSourceLonAndLatAtAlt2(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat,Double_t &sourceAltitude, Int_t maxLoopIter = -1, TVector3* sourcePosOptional = NULL) const;



  TVector3 getUnitVectorAlongThetaWavePhiWave(double thetaWave, double phiWave) const;
  
  int getSourceLonAndLatAtAlt3(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat,Double_t &sourceAltitude,
			       double* deltaAltIfNoIntersectection = NULL, bool returnBestPositionIfNoIntersection = false) const;
  /** 
   * New version of traceBackToContinent with same arguments as Cosmin's version.
   * However this one uses the new and lovely getSourceLonAndLatAtAlt3.
   * 
   * @param phiWave payload phi in radians (relative to adu5 aft-fore)
   * @param thetaWave radians, +ve theta is up, the silly EventCorrelator convention
   * @param lon pointer to the calculated longitude
   * @param lat pointer to the calculated latitude
   * @param alt pointer to the calculated altitud
   * @param theta_adjustment_required if non-null, will encode difference in theta required for intersection with the continent at the returned lon/lat/alt
   * 
   * @return same as getSourceLonAndLatAtAlt3
   */
  int traceBackToContinent3(Double_t phiWave, Double_t thetaWave, 
			    Double_t * lon, Double_t * lat, Double_t *alt, Double_t * theta_adjustment_required) const;


  
  

  /** 
   * Trace back to continent, based a bit on Abby's code. 
   * phiWave and thetaWave are in payload coordinates, and I believe should be in radians (they're just passed to getSourceLonAndLatAtDesiredAlt)
   *
   * This works by iterating altitude down until it hits the ground. If it never hits the ground, we try adjusting theta slightly, up to max_theta_adjustment, until it does. 
   * 
   * Returns 0 if never hits the ground, even with maximum adjustment
   *
   * Returns 1 if hits the ground with no adjustment
   *
   * Returns 2 if it hits the ground with adjustment
   *
   * If the pointers are passed, they are filled if we return 1 or 2, but not if we return 0. 
   *
   * A binary search with up to max_iter iterations is used to search for the smallest theta that hits the ground
   */ 
  int traceBackToContinent(Double_t phiWave, Double_t thetaWave, 
                           Double_t * lon, Double_t * lat, Double_t *alt, Double_t * theta_adjustment_required, 
                           Double_t max_theta_adjustment = TMath::Pi()/180, Int_t max_iter = 10) const; 




  /**
   * For a given source latitude, longitude and altitude calculates the payload centric azimuthal and elevation angles of the plane wave incident at the payload.
   * 
   * @param sourceLon The longitude of the source.
   * @param sourceLat The latitude of the source.
   * @param sourceAlt The altitude of the source.
   * @param phiWave Reference to a Double_t in which to store the azimuthal angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   * @param thetaWave Reference to a Double_t in which to store the elevation angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   */
  void getThetaAndPhiWave(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t &thetaWave, Double_t &phiWave);


  /**
   * const/thread-safe version of getThetaAndPhiWave, unlike that function the fSourcePos vector is not updated
   * 
   * @param sourceLon The longitude of the source.
   * @param sourceLat The latitude of the source.
   * @param sourceAlt The altitude of the source.
   * @param phiWave Reference to a Double_t in which to store the azimuthal angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   * @param thetaWave Reference to a Double_t in which to store the elevation angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   * @param sourcePos optional pointer to vector to store source position in
   */
  void getThetaAndPhiWave2(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t &thetaWave, Double_t &phiWave, TVector3* sourcePos = NULL) const;

  /** Cartesian version of above. Note that sourcePos will be modified so that it is rotated to balloon coords, so pass a copy if you don't want to change it */ 
  void getThetaAndPhiWaveCart(TVector3 * sourcePos, Double_t & thetaWave, Double_t & phiWave) const; 

  /** Find the direction in payload coordinates of a cartesian ray p0 - v0 t as t-> infinity or -infinity
   * This does it the dumbest possible way, there's probably a smart way to do it. 
   * */ 
  void getThetaAndPhiWaveOfRayAtInfinity(const TVector3 & p0, const TVector3 & v0 , Double_t & thetaWave, Double_t & phiWave,
                                          Bool_t plus_infinity = true, Double_t eps = 0.00001 * TMath::DegToRad(), Double_t step = 10e7, TVector3 * testPosition = 0) const; 

  /**
   * Returns the expected theta and phi expected from WAIS divide
   * 
   * @param phiWave Reference to a Double_t in which to store the azimuthal angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction of the ADU5 fore antenna)
   * @param thetaWave Reference to a Double_t in which to store the elevation angle of plane wave, in the convention of this function +ve theta is upwards.  
   */
  void getThetaAndPhiWaveWaisDivide(Double_t &thetaWave, Double_t &phiWave);


  
  /** 
   * Returns the expected theta and phi expected from LDB camp
   * 
   * @param phiWave Reference to a Double_t in which to store the azimuthal angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction of the ADU5 fore antenna)
   * @param thetaWave Reference to a Double_t in which to store the elevation angle of plane wave, in the convention of this function +ve theta is upwards.  
   */
  void getThetaAndPhiWaveLDB(Double_t &thetaWave, Double_t &phiWave);

  /** 
   * For a given base, calculates the theta angle from the base to the balloon - use as a horizon check
   * 
   * @param thetaWave Reference to a Double_t in which to store the elevation angle of plane wave from base to balloon  
   */
  void getThetaWaveAtBase(Double_t baseLon, Double_t baseLat, Double_t baseAlt, Double_t &thetaWave);


  
  /**
   * For a the Williams Field seavey antenna calculates the payload centric azimuthal and elevation angles of the plane wave incident at the payload.
   * 
   * @param phiWave Reference to a Double_t in which to store the azimuthal angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   * @param thetaWave Reference to a Double_t in which to store the elevation angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)  
   */
  void getThetaAndPhiWaveWillySeavey(Double_t &thetaWave, Double_t &phiWave);

  
  /*  For a the Williams Field borehole antenna calculates the payload centric azimuthal and elevation angles of the plane wave incident at the payload.
   * @param phiWave Reference to a Double_t in which to store the azimuthal angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   * @param thetaWave Reference to a Double_t in which to store the elevation angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   */
  void getThetaAndPhiWaveWillyBorehole(Double_t &thetaWave,Double_t &phiWave);

  
  /**
   * For a the Taylor Dome seavey antenna calculates the payload centric azimuthal and elevation angles of the plane wave incident at the payload.
   * 
   * @param phiWave Reference to a Double_t in which to store the azimuthal angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   * @param thetaWave Reference to a Double_t in which to store the elevation angle of plane wave (in payload centric coordinates with phi equals zero lying in the direction the ADU5 fore antenna)
   * */
  void getThetaAndPhiWaveTaylorDome(Double_t &thetaWave, Double_t &phiWave);


  
  /**
   * For a given source latitude, longitude and altitude calculates the plane wave crossing time difference between two antennas.
   * 
   * @param ant1 The first antenna.
   * @param ant2 The second antenna.
   * @param sourceLon The longitude of the source.
   * @param sourceLat The latitude of the source.
   * @param sourceAlt The altitude of the source.
   * @return The plane wave crossing time difference between the two antennas (t_1-t_2).
   */
  Double_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt); 


  
  /**
   * Calculates the plane wave crossing time difference between two antennas for a given wave phi and theta - useful for testing prior to flight when heading etc not known/needed
   * 
   * @param ant1 The first antenna.
   * @param ant2 The second antenna.
   * @return The plane wave crossing time difference between the two antennas (t_1-t_2).
   */
  Double_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave);

  /**
   * Calculates the plane wave crossing time difference between two antennas for a given wave cos and sin of phi and theta
   * useful for testing prior to flight when heading etc not known/needed, quicker than the above method if calculating deltaT for many phi and theta values
   * - no need to get cos and sin each time
   * 
   * @param ant1 The first antenna.
   * @param ant2 The second antenna.
   * @return The plane wave crossing time difference between the two antennas (t_1-t_2).
   */
  Double_t getDeltaTExpected(Int_t ant1,Int_t ant2,Double_t cosPhi,Double_t sinPhi,Double_t cosTheta,Double_t sinTheta);


  
  /**
   * Calculates the plane wave crossing time difference between two antennas for the Williams Field seavey.
   * 
   * @param ant1 The first antenna.
   * @param ant2 The second antenna.
   * @return The plane wave crossing time difference between the two antennas (t_1-t_2).
   */
  Double_t getDeltaTWillySeavey(Int_t ant1, Int_t ant2);


  
  /**
   * Calculates the plane wave crossing time difference between two antennas for the Williams Field borehole.
   * 
   * @param ant1 The first antenna.
   * @param ant2 The second antenna.
   * @return The plane wave crossing time difference between the two antennas (t_1-t_2).
   */
  Double_t getDeltaTWillyBorehole(Int_t ant1, Int_t ant2);


  
  /**
   * Calculates the plane wave crossing time difference between two antennas for the Taylor Dome borehole.
   * 
   * @param ant1 The first antenna.
   * @param ant2 The second antenna.
   * @return The plane wave crossing time difference between the two antennas (t_1-t_2).
   */
  Double_t getDeltaTTaylor(Int_t ant1, Int_t ant2);
  


  /** 
   * Payload centric with phi equals zero lying along the direction of the ADU5 fore antenna
   * 
   * @return Returns the last calculated azimuthal angle in radians.
   */
  Double_t getPhiWave() const {
    return fPhiWave;
  }

  /** 
   * Get last calculated elevation angle
   * @return the (payload centric) elevation angle last calculated.
   */
  Double_t getThetaWave() const {
    return fThetaWave;
  }

  /** 
   * Get the last calculated source longitude
   * @return the last calculated longitude
   */
  Double_t getSourceLongitude() const {
    return fSourceLongitude;
  }
  /** 
   * Get the last calculated source latitude
   * @return the last calculated source latitude
   */
  Double_t getSourceLatitude() const{
    return fSourceLatitude;
  }

  /** 
   * Get the last calculated source altitude.
   * @return the last calculated source altitude
   */  
  Double_t getSourceAltitude() const {
    return fSourceAltitude;
  }

  UInt_t getTaylorDomeTriggerTimeNs() const;								/// Gets the time of flight to Taylor Dome
  UInt_t getWaisDivideTriggerTimeNs() const;								/// Gets the time of flight to Wais Divide
  UInt_t getSipleTriggerTimeNs() const;									/// Gets the time of flight to Siple
  UInt_t getLDBTriggerTimeNs() const;									/// Gets the time of flight to LDB camp
  UInt_t getTriggerTimeNsFromSource(Double_t sourceLat, Double_t sourceLong, Double_t sourceAlt) const; /// Gets time of flight from any source
  Double_t getDistanceFromSource(Double_t sourceLat, Double_t sourceLong, Double_t sourceAlt) const;	/// Gets distance from any source in meters
  void setIncludeGroupDelay(Int_t flag){fIncludeGroupDelay=flag;}					/// Toggles the silly group delay correction on and off

  RampdemReader *fRampdemReader;
  AnitaGeomTool * fUPGeomTool; 

  Double_t getAngleBetweenPayloadAndSource(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt); //ACG additional 
  void getSunPosition(Double_t& phiDeg, Double_t& thetaDeg) const;
  Double_t getAzimuthOfSunRelativeToNorth() const;
  Double_t getAzimuthOfSun() const;
  Double_t getDifferencePointingToSun(Double_t phiAngle, Bool_t inputInDegrees=true) const;

  void getThetaAndPhiWaveHiCal(Double_t& thetaWave, Double_t& phiWave); /// Currently just does direct events...


  /** 
   * Returns a source elevation angle, assuming the image is a specular reflection of a source at infinity
   * no refraction is considered
   * added by S. Stafford 02/01/2017
   * @param plAlt payload altitude
   * @param imageEl elevation of the image in payload coords
   * @param imageAlt altitude at the reconstruction on continent 
   * 
   * @return Source elevation angle assuming specular reflection
   */
  static Double_t getReflectionAngle(Double_t plAlt, Double_t imageEl, Double_t imageAlt);

  /** 
   * Returns a source elevation angle, assuming the image is a specular reflection of a source at infinity
   * no refraction is considered
   * added by S. Stafford 02/01/2017
   * @param imageEl elevation of the image in payload coords
   * @param imageAlt altitude at the reconstruction on continent 
   * 
   * @return Source elevation angle assuming specular reflection
   */  
  Double_t getReflectionAngle(Double_t imageEl, Double_t imageAlt) const;


  /** 
   * Gets the astronomical coordinates corresponding to a wave coming from phiWave, thetaWave 
   * The sign convention for thetaWave is the same as elsewhere, where thetaWave is positive going down. 
   * Computed RA is in hours, dec, l and b in degrees. 
   * 
   * @param phiWave is the plane wave phi
   * @param thetaWave is the plane wave theta, with the anitaEventCorrelator sign convention
   * @param RA Right ascention in hours
   * @param dec declination in degrees
   * @param l (degrees)
   * @param b (degrees)
   * 
   * @return 
   */
  int astronomicalCoordinates(Double_t phiWave, Double_t thetaWave, Double_t * RA = 0, Double_t * dec = 0, Double_t * l = 0, Double_t * b = 0) const; 
  int fromRADec (Double_t RA,  Double_t dec, Double_t *phi, Double_t *theta) const; 
  

  /** 
   * Sets whether to print the debug messages
   * @param db true or false
   */
  inline void setDebug(bool db){fDebug = db;}

  /** 
   * Is the debug output on?
   * @return the value of fDebug
   */
  inline bool getDebug() const {return fDebug;}

  /** 
   * Set whether or not to use to interpolated version of RampdemReader::SurfaceAboveGeoid?
   * @param i true means use the bilinear interpolated version, false means use the plain old version
   */
  inline void setInterpSurfaceAboveGeoid(bool i){fInterpSurfaceAboveGeoid = i;} // *TOGGLE *GETTER=getInterpSurfaceAboveGeoid
  /** 
   * Is the flag set to use to bilinear interpolated version of RampdemReader::SurfaceAboveGeoid?
   * @return the value of fInterpSurfaceAboveGeoid
   */
  inline bool getInterpSurfaceAboveGeoid() const {return fInterpSurfaceAboveGeoid;}

  /** 
   * Sets the epsilon parameter at which the getSourceAtLonAndLatAtAlt thinks it's found a valid solution
   * 
   * @param closeEnough epsilon distance in metres (default is 1.0)
   */
  inline void setSurfaceCloseEnoughInter(double closeEnough = 1.0){fSurfaceCloseEnoughIter = closeEnough;}

  /** 
   * What is the currently set value for which getSourceAtLonAndLatAtAlt accepts agreement between source height and surfaceAboveGeoid?
   * @return the current value of fSurfaceCloseEnoughIter
   */
  inline double getSurfaceCloseEnoughInter() const {return fSurfaceCloseEnoughIter;}


  /** 
   * Set how many loop iterations to allow before giving up in getSourceLonAndLatAtAlt
   * 
   * @param n the number of loops to allow (default is 50, which should be enough in most instances)
   */
  inline void setMaxLoopIterations(Int_t n = 50){fMaxLoopIterations = n;}

  /** 
   * Get how many loop iterations to allow before giving up in getSourceLonAndLatAtAlt
   * @return the current value of fMaxLoopIterations 
   */
  inline int getMaxLoopIterations() const {return fMaxLoopIterations;}
  

  /** 
   * Get surface above geoid, depending on whether interpolation is requested.
   * @return surface elevation of the ice, above the geoid
   */
  inline double surfaceAboveGeoid(Double_t lon, Double_t lat) const {
    if(fInterpSurfaceAboveGeoid){
      return RampdemReader::BilinearInterpolatedSurfaceAboveGeoid(lon, lat);
    }
    else{
      return RampdemReader::SurfaceAboveGeoid(lon, lat);
    }
  }


  /**
   * Call this function if you update the payload lon/lat/alt after construction.
   * This function checks whether public longitude/latitude/altitude has changed
   * since the fBalloon cartesian stuff was calculated. Updates it if it has.
   */
  void updateCartesianBalloonInfo();

  // mutable std::vector<TGraph*> grTests;

 private:
  Int_t		fIncludeGroupDelay;		/// Include group delay in deltaTs? (default is no)
  TVector3	fSourcePos;			/// Private variable to hold the source location in cartesian coordinates.
  Double_t	fSourceLongitude;		/// The source longitude.
  Double_t	fSourceLatitude;		/// The source latitude.
  Double_t	fSourceAltitude;		/// The source altitude.
  Double_t	fThetaWave;			/// The elevation angle of the plane wave in payload centric coordinates.
  Double_t	fPhiWave;			/// The azimuthal angle of the plane wave in payload centric coordinates with phi equals zero lying along the direction of the ADU5 fore antenna.
  Double_t	fBalloonCoords[3];		/// The balloon position in cartesian coords
  TVector3	fBalloonPos;			/// The cartesian coords as a TVector3
  Double_t	fBalloonTheta;			/// The balloon theta
  Double_t	fBalloonPhi;			/// The balloon phi
  Double_t	fBalloonHeight;			/// The balloon height
  Float_t       fBalloonLonCache;		/// The public member longitude when the rest of the balloon coords were calculated
  Float_t       fBalloonLatCache;		/// The public member latitude when the rest of the balloon coords were calculated
  Float_t       fBalloonAltCache;		/// The public member altitude when the rest of the balloon coords were calculated
  Bool_t	fDebug;				/// Print lots of info in complicated methods, for debugging.
  Bool_t        fInterpSurfaceAboveGeoid;	/// Flag to use bilinear interpolation version of SurfaceAboveGeoid function (default = false)
  Double_t      fSurfaceCloseEnoughIter;	/// How close is close enough in the getSourceLonAndLatAtAlt iteration? (default = 1.0 metres)
  Int_t         fMaxLoopIterations;             /// How many iterations before giving up in getSourceLonAndLatAtAlt? (default = 50)
  ClassDef(UsefulAdu5Pat,0);			/// ROOT's magic macro.

  //optimisation stuff
  Double_t getDeltaTTaylorOpt(Int_t ant1, Int_t ant2, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi);
  Double_t getDeltaTExpectedOpt(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi);

  Double_t getDeltaTSeaveyOpt(Int_t ant1, Int_t ant2, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi);
  Double_t getDeltaTExpectedSeaveyOpt(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi);
  Double_t getGroupDelay(Double_t phiToAntBoresight);  
  

  /** 
   * Rotate the coordinate system so that theta/phi account for payload pitch and roll
   * This is a modularization of pre-existing code littered throughout this class
   * @param thetaWave updated after pitch/roll accounted for
   * @param phiWave updated after pitch/roll accounted for
   */
  void accountForPitchAndRollInPhiWaveThetaWave(Double_t& phiWave, Double_t& thetaWave) const;
  
};


#endif //USEFULADU5PAT_H
