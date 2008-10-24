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
#include "AnitaGeomTool.h"

//!  This is the position and attitude class that inherits from Adu5Pat and has some useful methods for determining the expected plane wave angles and inter-antenna crossing times for a given source.
/*!
  As UsefulAnitaEvent is to RawAnitaEvent, UsefulAdu5Pat is to Adu5Pat. Well not quite as useful but you get the picture.
*/
class UsefulAdu5Pat: public Adu5Pat
{

 public:
  UsefulAdu5Pat(); ///< Default constructor
  UsefulAdu5Pat(Adu5Pat *patPtr); ///< Assignment constructor
  ~UsefulAdu5Pat(); ///< Destructor

  //! For a given azimuthal and elevation angle of a plane wave (in payload coordinates) calculates the point on the Earth's surface that the source would come from.
  /*!
    \param phiWave Azimuthal angle of plane wave (in payload centric coordinates)
    \param thetaWave Elevation angle of plane wave (in payload centric coordinates)
    \param sourceLon Reference to a Double_t in which the longitude of the source will be stored.
    \param sourceLat Reference to a Double_t in which the latitude of the source will be stored.
    \return 1 on successful source localisation, o if the solution does not point to the ground.
  */
  int getSourceLonAndLatAltZero(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat);
  
//! For a given source latitude, longitude and altitude calculates the payload centric azimuthal and elevation angles of the plane wave incident at the payload.
  /*!
    
    \param sourceLon The longitude of the source.
    \param sourceLat The latitude of the source.
    \param sourceAlt The altitude of the source.
    \param phiWave Reference to a Fouble_t in which to store the azimuthal angle of plane wave (in payload centric coordinates)
    \param thetaWave Reference to a Fouble_t in which to store the elevation angle of plane wave (in payload centric coordinates)
  */
  void getThetaAndPhiWave(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t &thetaWave, Double_t &phiWave);
  //! For a the Williams Field seavey antenna calculates the payload centric azimuthal and elevation angles of the plane wave incident at the payload.
  /*!    
    \param phiWave Reference to a Fouble_t in which to store the azimuthal angle of plane wave (in payload centric coordinates)
    \param thetaWave Reference to a Fouble_t in which to store the elevation angle of plane wave (in payload centric coordinates)
  */
  void getThetaAndPhiWaveWillySeavey(Double_t &thetaWave, Double_t &phiWave);
  //! For a given source latitude, longitude and altitude calculates the plane wave crossing time difference between two antennas.
  /*!
    \param ant1 The first antenna.
    \param ant2 The second antenna.
    \param sourceLon The longitude of the source.
    \param sourceLat The latitude of the source.
    \param sourceAlt The altitude of the source.
    \return The plane wave crossing time difference between the two antennas (t_1-t_2).
  */
  Double_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt); 
   //! Calculates the plane wave crossing time difference between two antennas for the Williams Field seavey.
  /*!
    \param ant1 The first antenna.
    \param ant2 The second antenna.
    \return The plane wave crossing time difference between the two antennas (t_1-t_2).
  */
  Double_t getDeltaTWillySeavey(Int_t ant1, Int_t ant2); 
  //! Calculates the plane wave crossing time difference between two antennas for the Williams Field borehole.
  /*!
    \param ant1 The first antenna.
    \param ant2 The second antenna.
    \return The plane wave crossing time difference between the two antennas (t_1-t_2).
  */
  Double_t getDeltaTWillyBorehole(Int_t ant1, Int_t ant2);

  
  Double_t getPhiWave() 
    { return fPhiWave;} ///< Returns the (payload centric) azimuthal angle last calculated.
  Double_t getThetaWave() 
     { return fThetaWave;} ///< Returns the (payload centric) elevation angle last calculated.
  Double_t getSourceLongitude() 
     { return fSourceLongitude;} ///< Returns the last calculated source longitude.
  Double_t getSourceLatitude() 
     { return fSourceLatitude;} ///< Returns the last calculated source latitude.
  Double_t getSourceAltitude() 
     { return fSourceAltitude;} ///< Returns the last calculated source altitude.
 private:
  TVector3 fSourcePos; ///< Private variable to hold the source location in cartesian coordinates.
  Double_t fSourceLongitude; ///< The source longitude.
  Double_t fSourceLatitude; ///< The source latitude.
  Double_t fSourceAltitude; ///< The source altitude.
  Double_t fThetaWave; ///< The elevation angle of the plane wave in payload centric coordinates.
  Double_t fPhiWave; ///< The azimuthal angle of the plane wave in payload centric coordinates.

  ClassDef(UsefulAdu5Pat,1); ///< ROOT's magic macro.
};


#endif //USEFULADU5PAT_H
