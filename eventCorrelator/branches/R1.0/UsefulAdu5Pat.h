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

class UsefulAdu5Pat: public Adu5Pat
{

 public:
  UsefulAdu5Pat();
  UsefulAdu5Pat(Adu5Pat *patPtr);
  ~UsefulAdu5Pat();


  int getSourceLonAndLatAltZero(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat);
  void getThetaAndPhiWave(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t &thetaWave, Double_t &phiWave);
  void getThetaAndPhiWaveWillySeavey(Double_t &thetaWave, Double_t &phiWave);
  Double_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt); 
  Double_t getDeltaTWillySeavey(Int_t ant1, Int_t ant2);
  Double_t getDeltaTWillyBorehole(Int_t ant1, Int_t ant2);
  Double_t getPhiWave() 
     { return fPhiWave;}
  Double_t getThetaWave() 
     { return fThetaWave;}
  Double_t getSourceLongitude() 
     { return fSourceLongitude;}
  Double_t getSourceLatitude() 
     { return fSourceLatitude;}
  Double_t getSourceAltitude() 
     { return fSourceAltitude;}
 private:
  TVector3 fSourcePos;
  Double_t fSourceLongitude;
  Double_t fSourceLatitude;
  Double_t fSourceAltitude;
  Double_t fThetaWave;
  Double_t fPhiWave;

  ClassDef(UsefulAdu5Pat,1);
};


#endif //USEFULADU5PAT_H
