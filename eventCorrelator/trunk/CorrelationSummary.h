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
  
  //!  There are 35 correlations formed from a set of 13 antennas (5-phi sectors).
   /*!
     For clarity in the list below the antennas are described using the following names, for cases where there is a nadir antenna in the centre phi sector we have:
    - LLT LT CT RT RRT
    - LLB LB CB RB RRB
    - LLN X  CN X  RRN
    In this scheme the 35 correlations are:
    -  "[0]" LT - LB
    -  "[1]" CT - CB
    -  "[2]" RT - RB
    -  "[3]" LT - CT
    -  "[4]" CT - RT
    -  "[5]" LB - CB
    -  "[6]" CB - RB
    -  "[7]" LT - CB
    -  "[8]" RT - CB
    -  "[9]" CT - LB
    -  "[10]" CT - RB
    -  "[11]" LLT - CT
    -  "[12]" CT - RRT
    -  "[13]" LLB - CB
    -  "[14]" CB - RRB
    -  "[15]" LLT - LT
    -  "[16]" RT - RRT
    -  "[17]" LLB - LB
    -  "[18]" RB - RRT 
    -  "[19]" CN - LLN
    -  "[20]" CN - RRN
    -  "[21]" LLN - LLB
    -  "[22]" LLN - LB
    -  "[23]" CN - CB
    -  "[24]" CN - LB
    -  "[25]" CN - RB
    -  "[26]" RRN - RB
    -  "[27]" RRN - RRB
    -  "[28]" LLN - LLT
    -  "[29]" LLN - LT
    -  "[30]" CN - CT
    -  "[31]" CN - LT
    -  "[32]" CN - RT
    -  "[33]" RRN - RT
    -  "[34]" RRN - RRT
    

    For cases where there isn't a nadir antenna in the centre phi sector, then we have:
     - LLT LT CT RT RRT
    - LLB LB CB RB RRB
    -  X  LN X  RN  X
    In this scheme the 35 correlations are (RN-LT is missing):
    -  "[0]" LT - LB
    -  "[1]" CT - CB
    -  "[2]" RT - RB
    -  "[3]" LT - CT
    -  "[4]" CT - RT
    -  "[5]" LB - CB
    -  "[6]" CB - RB
    -  "[7]" LT - CB
    -  "[8]" RT - CB
    -  "[9]" CT - LB
    -  "[10]" CT - RB
    -  "[11]" LLT - CT
    -  "[12]" CT - RRT
    -  "[13]" LLB - CB
    -  "[14]" CB - RRB
    -  "[15]" LLT - LT
    -  "[16]" RT - RRT
    -  "[17]" LLB - LB
    -  "[18]" RB - RRT 
    -  "[19]" LN - RN
    -  "[20]" LN - LLB
    -  "[21]" LN - LB
    -  "[22]" LN - CB
    -  "[23]" LN - RB
    -  "[24]" RN - RRB
    -  "[25]" RN - RB
    -  "[26]" RN - CB
    -  "[27]" RN - LB
    -  "[28]" LN - LLT
    -  "[29]" LN - LT
    -  "[30]" LN - CT
    -  "[31]" LN - RT
    -  "[32]" RN - RRT
    -  "[33]" RN - RT
    -  "[34]" RN - CT
   */
  int firstAnt[35]; 
  int secondAnt[35]; ///< The index of the second antenna in the 35 possible pairs (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).
  double maxCorVals[35]; ///< The maximumn correlation value for each of the 35 possible correlations (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).
  double maxCorTimes[35]; ///< The time of the maximum correlation value for each of the 35 possible correlations (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).
  double rmsCorVals[35]; ///< The rms correlation value for each of the 35 possible correlations (3 up-down, 4 left-right, 4 diagonal, 4 outside-neighbour).




  double secondCorVals[35][2]; ///< The peak of the next highest correlation values (tore both left and right vals).
  double secondCorTimes[35][2]; ///< The time of the next highest correlation values (tore both left and right vals).

  //Time Thingies
  //  double deltaT[35]; is maxCorTimes
  double deltaTErr[35]; ///< The error on each of the 35 deltaTs (no idea right now what this means).
  
  //Fit Results
  double phiWave; ///< The azimuthal angle of the plane wave (in payload centric coordinates).
  double thetaWave; ///< The elevation angle of the plane wave (in payload centric coordinates).
  double phiWaveErr; ///< The error on the azimuthal angle of the plane wave (in payload centric coordinates).
  double thetaWaveErr; ///< The error on the elevation angle of the plane wave (in payload centric coordinates).
  double chiSq; ///< The chi-squared of the fit.
  int ndf; ///< The number of degrees of freedom -- no frigging idea I just make it up

  //Antenna postion variables for use in fit
  Double_t fAntPhi[35][2]; ///< A lookup table for antenna postions.
  Double_t fAntR[35][2]; ///< A lookup table for antenna postions.
  Double_t fAntZ[35][2]; ///< A lookup table for antenna postions.
  
  


  ClassDef(CorrelationSummary,4); ///< One of ROOT's magic macros.
};


#endif //CORRELATIONSUMMARY_H
