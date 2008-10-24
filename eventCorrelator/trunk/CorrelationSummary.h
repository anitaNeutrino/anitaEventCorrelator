///////////////////////////////////////////////////////////////////////////////
/////  CorrelationSummary.h        ANITA Correlation Summary              /////
/////                                                                     /////
/////  Description:                                                       /////
/////     A simple class for storing correlation summary info for an event/////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                            /////
///////////////////////////////////////////////////////////////////////////////

#ifndef CORRELATIONSUMMARY_H
#define CORRELATIONSUMMARY_H

#include "TObject.h"
#include "TGraph.h"


//!  This is a poorly thought out class that was meant to be a summary of peaks of the correlations of soem antennas.
/*!
  This is a poorly thought outclass that was meant to be a summary of peaks of the correlations of soem antennas.
*/
class CorrelationSummary: public TObject
{


 public:
   CorrelationSummary(); ///< Default constructor
   ~CorrelationSummary(); ///< Destructor
   
  //! Assignment constructor.
  /*!
    \param teventNumber The event number.
    \param tcentreAnt The centre antenna used in the correlation.
    \param sixAnts An array of the six antennas used.
    \param deltaT The sampling period used in the interpolation (or zero)    
  */
  CorrelationSummary(int teventNumber, int tcentreAnt, int sixAnts[6], double deltaT=0);
  void fillErrorsAndFit(); ///< The worker function that actually does the fitting.

  //! Tests a given plane wave hypothesis using a either six or ten antennas.
  /*!
    \param tPhiWave The (payload centric) azimuthal angle of the plane wave.
    \param tThetaWave The (payload centric) elevation angle of the plane wave.
    \param numAnts The number of antennas to use in the test.
    \return The chi-squared of this plane wave hypothesis.
  */
  Double_t getChiSquared(Double_t tPhiWave, Double_t tThetaWave, Int_t numAnts);
  
  //! For a given plane wave hypothesis returns the expected time difference between one of the pairs of antennas.
  /*!
    \param tPhiWave The (payload centric) azimuthal angle of the plane wave.
    \param tThetaWave The (payload centric) elevation angle of the plane wave.
    \param pairInd The index of the pair (the numbering system is the 3 top bottom pairs, 4 left-right pairs, 4 diagonal pairs, plus the 4 'neighbour' + the 4 (next phi to neighbour).
    \return The expected deltaT for this plane wave and this pair of antennas.
  */
  Double_t getDeltaTExpected(Double_t tPhiWave, Double_t tThetaWave, Int_t pairInd);
  
  void setFitResults(Double_t tPhi, Double_t tTheta, Double_t tPhiErr, Double_t tThetaErr, Double_t tChiSq); ///< Sets the result of an external fit.
 

  //Simple Event Characteristics
  int eventNumber; ///< The eventNumber.
  int centreAntenna; ///< The number of the centre antenna.
  int sixAnts[6]; ///< The numbers of the six central antennas.
  int nextFourAnts[4]; ///< The numbers of the four outside antennas.
  double deltaT; ///< The sampling period used.


  //Correlation Thingies
  int firstAnt[19]; ///< The index of the first antenna in the 19 possible pairs (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).
  int secondAnt[19]; ///< The index of the second antenna in the 19 possible pairs (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).
  double maxCorVals[19]; ///< The maximumn correlation value for each of the 19 possible correlations (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).
  double maxCorTimes[19]; ///< The time of the maximum correlation value for each of the 19 possible correlations (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).
  double rmsCorVals[19]; ///< The rms correlation value for each of the 19 possible correlations (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).
s 


  double secondCorVals[19][2]; ///< The peak of the next highest correlation values (tore both left and right vals).
  double secondCorTimes[19][2]; ///< The time of the next highest correlation values (tore both left and right vals).

  //Time Thingies
  //  double deltaT[19]; is maxCorTimes
  double deltaTErr[19]; ///< The error on each of the 19 deltaTs (no idea right now what this means).
  
  //Fit Results
  double phiWave; ///< The azimuthal angle of the plane wave (in payload centric coordinates).
  double thetaWave; ///< The elevation angle of the plane wave (in payload centric coordinates).
  double phiWaveErr; ///< The error on the azimuthal angle of the plane wave (in payload centric coordinates).
  double thetaWaveErr; ///< The error on the elevation angle of the plane wave (in payload centric coordinates).
  double chiSq; ///< The chi-squared of the fit.
  int ndf; ///< The number of degrees of freedom -- no frigging idea I just make it up

  //Antenna postion variables for use in fit
  Double_t fAntPhi[19][2]; ///< A lookup table for antenna postions.
  Double_t fAntR[19][2]; ///< A lookup table for antenna postions.
  Double_t fAntZ[19][2]; ///< A lookup table for antenna postions.
  
  


  ClassDef(CorrelationSummary,3); ///< One of ROOT's magic macros.
};


#endif //CORRELATIONSUMMARY_H
