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
  
  


  ClassDef(CorrelationSummary,2);
};


#endif //CORRELATIONSUMMARY_H
