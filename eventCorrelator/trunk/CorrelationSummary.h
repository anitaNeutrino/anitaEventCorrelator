//////////////////////////////////////////////////////////////////////////////
/////  CorrelationSummary.h        ANITA Correlation Summary              /////
/////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for storing correlation summary info for an event/////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#ifndef CORRELATIONSUMMARY_H
#define CORRELATIONSUMMARY_H

#include "TObject.h"
#include "TGraph.h"



class CorrelationSummary: public TObject
{


 public:
  CorrelationSummary(int teventNumber, int tmaxAnt, int sixAnts[6]);
  void fillErrorsAndFit();
  


  //Simple Event Characteristics
  int eventNumber;
  int maxAntenna;
  int sixAnts[6];


  //Correlation Thingies
  int firstAnt[11];
  int secondAnt[11];
  double maxCorVals[11]; //For the 3 top bottom pairs, 4 left-right pairs, 4 diagonal pairs 
  double maxCorTimes[11];
  double rmsCorVals[11]; //For the 3 top bottom pairs, 4 left-right pairs, 4 diagonal pairs 


  double secondCorVals[11][2]; //Store both left and right vals
  double secondCorTimes[11][2]; 

  //Time Thingies
  //  double deltaT[11]; is maxCorTimes
  double deltaTErr[11];
  
  //Fit Results
  double phiWave;
  double thetaWave;
  
  


  ClassDef(CorrelationSummary,1);
};


#endif //CORRELATIONSUMMARY_H
