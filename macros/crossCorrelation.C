//////////////////////////////////
///// CrossCorrelation.C
///// A function that does some cross correlations and
///// interferometric map calculations ... basics of this have been
///// stolen from Ryan
///// Author: Matthew Mottram
//////////////////////////////////

#include <iostream>
#include <fstream>

//ANITA Includes
#include "PrettyAnitaEvent.h"
#include "UsefulAdu5Pat.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "RawAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "FFTtools.h"
#include "FFTWComplex.h"


//ROOT Includes
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"

#define NUM_BINS_THETA 90
#define NUM_BINS_PHI 90
#define PI 3.14159265

//runs from palestine to use: 3847 - 3853, 3793, 3797, 3798

TGraph *getCorrelation(TGraph *gr1,TGraph *gr2);
void getTriggeredPhi(RawAnitaHeader *hdPtr,int triggeredPhi[16]);
void crossCorrelate(RawAnitaEvent *evPtr,RawAnitaHeader *hdPtr,Adu5Pat *patPtr);
void startCorrelation(int run,int entry);
void setupCosSinArray(double thetaArray[NUM_BINS_THETA],double phiArray[NUM_BINS_PHI],double cosThetaArray[NUM_BINS_THETA],double sinThetaArray[NUM_BINS_THETA],double cosPhiArray[NUM_BINS_PHI],double sinPhiArray[NUM_BINS_PHI]);

void startCorrelation(int run,int entry){

  char headName[FILENAME_MAX];
  char eventName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  RawAnitaEvent *evPtr = 0;
  RawAnitaHeader *hdPtr = 0;
  Adu5Pat *patPtr =0;

  sprintf(headName,"/data/anita/palestine08/root/run%d/headFile%d.root",run,run);
  sprintf(eventName,"/data/anita/palestine08/root/run%d/eventFile%d.root",run,run);
  sprintf(gpsName,"/data/anita/palestine08/root/run%d/gpsFile%d.root",run,run);

  TFile *eventFile = new TFile(eventName);
  TFile *headFile = new TFile(headName);
  TFile *fpGps = new TFile(gpsName);

  TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
  TTree *eventTree = (TTree*)eventFile->Get("eventTree");
  TTree *headTree = (TTree*)headFile->Get("headTree");

  eventTree->SetBranchAddress("event",&evPtr);
  headTree->SetBranchAddress("header",&hdPtr);
  adu5PatTree->SetBranchAddress("pat",&patPtr);

  eventTree->GetEntry(entry);
  headTree->GetEntry(entry);
  adu5PatTree->GetEntry(entry);

  while(hdPtr->triggerTimeNs<649.5e6 || hdPtr->triggerTimeNs>650.5e6){
    std::cout << "opened entry " << entry << " (event " << hdPtr->eventNumber << " with triggerTimeNs " << hdPtr->triggerTimeNs << std::endl;
    entry++;
    eventTree->GetEntry(entry);
    headTree->GetEntry(entry);
    adu5PatTree->GetEntry(entry);
  }
 std::cout << "opened entry " << entry << " (event " << hdPtr->eventNumber << " with triggerTimeNs " << hdPtr->triggerTimeNs << std::endl;
    entry++;
    crossCorrelate(evPtr,hdPtr,patPtr);
 
}


TGraph *getCorrelation(TGraph *gr1,TGraph *gr2){
  return FFTtools::getCorrelationGraph(gr1,gr2);
}


void getTriggeredPhi(RawAnitaHeader *hdPtr,int triggeredPhi[16]){
  for(int phi=0;phi<16;phi++){
    if(hdPtr->l3TrigPattern & (1 << phi)) triggeredPhi[phi]=1;
    else triggeredPhi[phi]=0;
  }
}


void getTriggeredAnt(int triggeredPhi[16],int triggeredAnt[32]){
  int phi;
  for(int ant=0;ant<32;ant++){
    phi=AnitaGeomTool::getPhiFromAnt(ant);
    if(triggeredPhi[phi]) triggeredAnt[ant]=1;
    else triggeredAnt[ant]=0;
    //std::cout << "ant " << ant << " phi " << phi << " triggered? " << triggeredAnt[ant] << " " << triggeredPhi[phi] << std::endl;
  }
}


void setupCosSinArray(double thetaArray[NUM_BINS_THETA],double phiArray[NUM_BINS_PHI],double cosThetaArray[NUM_BINS_THETA],double sinThetaArray[NUM_BINS_THETA],double cosPhiArray[NUM_BINS_PHI],double sinPhiArray[NUM_BINS_PHI]){
  for(int i=0;i<NUM_BINS_PHI;i++){
    phiArray[i] = (i+0.5) * 2*PI/NUM_BINS_PHI - PI;
    cosPhiArray[i] = cos(phiArray[i]);
    sinPhiArray[i] = sin(phiArray[i]);
  }
  for(int i=0;i<NUM_BINS_THETA;i++){
    thetaArray[i] = (i+0.5) * PI/NUM_BINS_THETA - PI/2;
    cosThetaArray[i] = cos(thetaArray[i]);
    sinThetaArray[i] = sin(thetaArray[i]);
  }
}


void crossCorrelate(RawAnitaEvent *evPtr,RawAnitaHeader *hdPtr,Adu5Pat *patPtr){

  gStyle->SetPalette(1);

  int triggeredPhi[16];
  int triggeredAnt[32];

  double thetaArray[NUM_BINS_THETA];
  double phiArray[NUM_BINS_PHI];
  double cosThetaArray[NUM_BINS_THETA];
  double sinThetaArray[NUM_BINS_THETA];
  double cosPhiArray[NUM_BINS_PHI];
  double sinPhiArray[NUM_BINS_PHI];

  PrettyAnitaEvent realEvent(evPtr,WaveCalType::kVTFullAGCrossCorClock,hdPtr);
  UsefulAdu5Pat usefulPat(patPtr);

  getTriggeredPhi(hdPtr,triggeredPhi);
  getTriggeredAnt(triggeredPhi,triggeredAnt);
  setupCosSinArray(thetaArray,phiArray,cosThetaArray,sinThetaArray,cosPhiArray,sinPhiArray);

  TGraph *grInt[32]={0};
  TGraph *grCorr[528]={0};
  double firstCorrTime[528];
  double dummyY;
  Double_t deltaT=1/(2.6*8);

  for(int ant=0;ant<32;ant++){
    if(triggeredAnt[ant]){
      grInt[ant]=realEvent.getGraph(AnitaGeomTool::getChanIndexFromAntPol(ant,AnitaPol::kVertical));
    }
  }
  int arrayRef=0;

  for(int ant1=0;ant1<32;ant1++){
    for(int ant2=ant1;ant2<32;ant2++){

      if(triggeredAnt[ant1] && triggeredAnt[ant2] && ant1!=ant2){
	grCorr[arrayRef]=FFTtools::getInterpolatedCorrelationGraph(grInt[ant1],grInt[ant2],deltaT);
	grCorr[arrayRef]->GetPoint(0,firstCorrTime[arrayRef],dummyY);
	arrayRef++;
      }
      else{
	if(ant1!=ant2) arrayRef++;
      }

    }
  }

  double deltaTarray[NUM_BINS_PHI][NUM_BINS_THETA];
  double correlationArray[NUM_BINS_PHI][NUM_BINS_THETA];

  for(int phi=0;phi<NUM_BINS_PHI;phi++){
    for(int theta=0;theta<NUM_BINS_THETA;theta++){
      deltaTarray[phi][theta]=0.;
      correlationArray[phi][theta]=0.;
    }
  }

  int getPoint;
  double xVal1;
  double xVal2;
  double weight1;
  double weight2;
  double pointVal1;
  double pointVal2;

  arrayRef=0;
  for(int ant1=0;ant1<32;ant1++){
    for(int ant2=ant1;ant2<32;ant2++){

      if(triggeredAnt[ant1] && triggeredAnt[ant2] && ant1!=ant2){// && ant1==2 && ant2==21){

	for(int phi=0;phi<NUM_BINS_PHI;phi++){
	  for(int theta=0;theta<NUM_BINS_THETA;theta++){

	    deltaTarray[phi][theta] = usefulPat.getDeltaTExpected(ant1,ant2,cosPhiArray[phi],sinPhiArray[phi],cosThetaArray[theta],sinThetaArray[theta]);
	    //deltaTarray[phi][theta] = usefulPat.getDeltaTExpected(ant1,ant2,phiArray[phi],thetaArray[theta]);

	    getPoint=static_cast<int>((deltaTarray[phi][theta]-firstCorrTime[arrayRef])*1/deltaT);

	    grCorr[arrayRef]->GetPoint(getPoint,xVal1,pointVal1);
	    grCorr[arrayRef]->GetPoint(getPoint+1,xVal2,pointVal2);

	    weight1 = 1 - fabs(deltaTarray[phi][theta]-xVal1)/(deltaT);
	    weight2 = 1 - fabs(deltaTarray[phi][theta]-xVal2)/(deltaT);

	    correlationArray[phi][theta]+=(pointVal1*weight1+pointVal2*weight2);

	  }//theta
	}//phi

      }//if

      if(ant1!=ant2) arrayRef++;

    }//ant2
  }//ant1

  char histName[FILENAME_MAX];
  sprintf(histName,"sumCrossCorrs");
  TH2D *sumCrossCorrs = new TH2D(histName,histName,NUM_BINS_PHI,-PI*180/PI,PI*180/PI,NUM_BINS_THETA,-PI/2*180/PI,PI/2*180/PI);
  
  for(int phi=0;phi<NUM_BINS_PHI;phi++){
    for(int theta=0;theta<NUM_BINS_THETA;theta++){
      sumCrossCorrs->Fill(phiArray[phi]*180/PI,thetaArray[theta]*180/PI,correlationArray[phi][theta]);
    }
  }

  sprintf(histName,"crossCorrCan");
  TCanvas *crossCorrCan = new TCanvas(histName,histName,800,400);
  sumCrossCorrs->Draw("aitoff");

}
