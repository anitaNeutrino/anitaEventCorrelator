//////////////////////////////////////////////////////////////////////////////
/////  CorrelationSummaryAnita3.cxx        ANITA event reading class     /////
/////  Description:                                                      /////
/////     A simple class for plotting event stuff like waveforms and     /////
/////     correlations adapted to the ANITA-3 geometry                   /////
/////  Author: Linda Cremonesi (l.cremonesi@ucl.ac.uk)                   /////
/////         based on CorrelarationSummary.cxx written by               ///// 
/////         Ryan Nichol (rjn@hep.ucl.ac.uk)                            /////
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//ANITA Includes
#include "CorrelationSummaryAnita3.h"
#include "AnitaGeomTool.h"

//ROOT Includes
#include "TMath.h"
#include "TStyle.h"
#include "TMinuit.h"


ClassImp(CorrelationSummaryAnita3);


//Global ish things.
AnitaGeomTool *fCSGeomToolAnita3=0;


CorrelationSummaryAnita3::CorrelationSummaryAnita3( )
{

   if(!fCSGeomToolAnita3)
      fCSGeomToolAnita3=AnitaGeomTool::Instance();
}

CorrelationSummaryAnita3::~CorrelationSummaryAnita3( )
{

}
CorrelationSummaryAnita3::CorrelationSummaryAnita3( int teventNumber, int tcentreAnt, int tnineAnts[], double dt)
   : eventNumber(teventNumber),centreAntenna(tcentreAnt),deltaT(dt)
{
  for(int i=0;i<9;i++) 
    nineAnts[i]=tnineAnts[i];

   if(!fCSGeomToolAnita3)
      fCSGeomToolAnita3=AnitaGeomTool::Instance();
}


void CorrelationSummaryAnita3::fillErrorsAndFit()
{
   //At some stage this should probably do something, otherwise it won't be much cop.

   //For now we'll go for the zeroth order solution to errors that is will assign them all to be the same
   //and just for teh sake of doing something we'll make it half a time bin (0.5/2.6) in ns.
   for(int i=0;i< NUM_CORRELATIONS_ANITA3;i++) {
      deltaTErr[i]=0.5/2.6; // in ns
   }

  //Now fill the antenna postions (these might become member variables)
   for(int i=0;i< NUM_CORRELATIONS_ANITA3;i++) {
      fAntPhi[i][0]=fCSGeomToolAnita3->getAntPhiPositionRelToAftFore(firstAnt[i]);
      fAntPhi[i][1]=fCSGeomToolAnita3->getAntPhiPositionRelToAftFore(secondAnt[i]);
      fAntR[i][0]=fCSGeomToolAnita3->getAntR(firstAnt[i]);
      fAntR[i][1]=fCSGeomToolAnita3->getAntR(secondAnt[i]);
      fAntZ[i][0]=fCSGeomToolAnita3->getAntZ(firstAnt[i]);
      fAntZ[i][1]=fCSGeomToolAnita3->getAntZ(secondAnt[i]);
   }
   

}

Double_t CorrelationSummaryAnita3::getChiSquared(Double_t tPhiWave, Double_t tThetaWave, Int_t numAnts)
{
   Double_t chiSq=0;
   for(int i=0;i<numAnts;i++) {
      Double_t dtExpect=getDeltaTExpected(tPhiWave,tThetaWave,i);
      chiSq+=(maxCorTimes[i]-dtExpect)*(maxCorTimes[i]-dtExpect)/(deltaTErr[i]*deltaTErr[i]);
      //      std::cout << i << "\t" << chiSq << "\t" << maxCorTimes[i] 
      //		<< "\t" << dtExpect << "\t"
      //		<<  deltaTErr[i] << "\n";
   }
   return chiSq;
}

Double_t CorrelationSummaryAnita3::getDeltaTExpected(Double_t tPhiWave, Double_t tThetaWave, Int_t pairInd)
{

   Double_t tanThetaW=TMath::Tan(tThetaWave);
   Double_t part1=fAntZ[pairInd][0]*tanThetaW - fAntR[pairInd][0] * TMath::Cos(tPhiWave-fAntPhi[pairInd][0]);
   Double_t part2=fAntZ[pairInd][1]*tanThetaW - fAntR[pairInd][1] * TMath::Cos(tPhiWave-fAntPhi[pairInd][1]);
   return  1e9*((TMath::Cos(tThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}

void CorrelationSummaryAnita3::setFitResults(Double_t tPhi, Double_t tTheta, Double_t tPhiErr, Double_t tThetaErr, Double_t tChiSq) {
   phiWave=tPhi;
   thetaWave=tTheta;
   phiWaveErr=tPhiErr;
   thetaWaveErr=tThetaErr;
   chiSq=tChiSq;
}
