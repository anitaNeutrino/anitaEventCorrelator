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



class CorrelationSummary: public TObject
{


 public:
   CorrelationSummary();
   ~CorrelationSummary();
  CorrelationSummary(int teventNumber, int tcentreAnt, int sixAnts[6], double deltaT=0);
  void fillErrorsAndFit();
  Double_t getChiSquared(Double_t tPhiWave, Double_t tThetaWave, Int_t numAnts);
  Double_t getDeltaTExpected(Double_t tPhiWave, Double_t tThetaWave, Int_t pairInd);
  void setFitResults(Double_t tPhi, Double_t tTheta, Double_t tPhiErr, Double_t tThetaErr, Double_t tChiSq);
 

  //Simple Event Characteristics
  int eventNumber;
  int centreAntenna;
  int sixAnts[6];
  int nextFourAnts[4];
  double deltaT;


  //Correlation Thingies
  int firstAnt[19];
  int secondAnt[19];
  double maxCorVals[19]; //For the 3 top bottom pairs, 4 left-right pairs, 4 diagonal pairs, plus the 4 'neighbour' + the 4 (next phi to neighbour)
  double maxCorTimes[19];
  double rmsCorVals[19]; //For the 3 top bottom pairs, 4 left-right pairs, 4 diagonal pairs 


  double secondCorVals[19][2]; //Store both left and right vals
  double secondCorTimes[19][2]; 

  //Time Thingies
  //  double deltaT[19]; is maxCorTimes
  double deltaTErr[19];
  
  //Fit Results
  double phiWave;
  double thetaWave;
  double phiWaveErr;
  double thetaWaveErr;
  double chiSq;
  int ndf; //no frigging idea

  //Antenna postion variables for use in fit
  Double_t fAntPhi[19][2];
  Double_t fAntR[19][2];
  Double_t fAntZ[19][2];
  
  


  ClassDef(CorrelationSummary,3);
};


#endif //CORRELATIONSUMMARY_H
