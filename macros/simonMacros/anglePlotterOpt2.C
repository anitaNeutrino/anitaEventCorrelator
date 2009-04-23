#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "UsefulAdu5Pat.h"
#include "CorrelationSummary.h"
#include "AnitaGeomTool.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>

double deltaTArrayMod[32];

void setArray();
void setArray(){
  for(int i = 0; i < 32; i++){
    deltaTArrayMod[i] = 0;
  }
}


void setValue(int x, double t);
void setValue(int x, double t){
  
    deltaTArrayMod[x] = t;
  
}



double anglePlotterOpt2(double *par);  
double anglePlotterOpt2(double *par){
  
  double theReturn = 99;
  double sumMean = 0;
  double sumMean2 = 0;

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corrName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  char baseDir[FILENAME_MAX];
  char *corTreeDir = "/home/swb/ANITA/eventCorr/outFiles";

  int run = 18;

  //sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
  sprintf(baseDir,"/unix/anita3/flight0809/root/");
  sprintf(eventName,"%s/run%d/eventFile%d.root",baseDir,run,run);
  sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
  sprintf(gpsName,"%s/run%d/gpsFile%d.root",baseDir,run,run);
  sprintf(corrName,"%s/corRun%d.root",corTreeDir,run);

  // TCanvas *canSurf= new TCanvas("canSurf","canSurf");
  // canSurf->Divide(1,2);

  vector<vector<double> > phiAngle;
  vector<vector<double> > deltaTVec;
  vector<vector<int> > firstAntVec;
  vector<vector<int> > secondAntVec;

  vector<double> temp;
  vector<int> temp2;
  temp.push_back(0);
  temp2.push_back(0);


  for(int i = 0; i < 32; i++){
    phiAngle.push_back(temp);
    deltaTVec.push_back(temp);
    firstAntVec.push_back(temp2);
    secondAntVec.push_back(temp2);
  }

 for(int i = 0; i < 16; i++){
   // cout << par[i] << endl;
 }

   
  RawAnitaEvent *event = 0;
  PrettyAnitaHk *hk = 0;
   
  RawAnitaHeader *header =0;
  Adu5Pat *pat =0;
  CorrelationSummary *corSum =0;
   
   
  TFile *fpHead = TFile::Open(headerName);
  TTree *headTree = (TTree*) fpHead->Get("headTree");
  headTree->SetBranchAddress("header",&header);
  headTree->BuildIndex("eventNumber");
   
  TFile *fpGps = TFile::Open(gpsName);
  TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
  adu5PatTree->BuildIndex("realTime");
  adu5PatTree->SetBranchAddress("pat",&pat);
   
  Int_t labChip;
  TFile *fpCor = new TFile(corrName);
  TTree *corTree = (TTree*) fpCor->Get("corTree");
  corTree->SetBranchAddress("cor",&corSum);
  corTree->SetBranchAddress("labChip",&labChip);

  Long64_t numEntries=corTree->GetEntries();
  int counter=0;

  Long64_t entry=0;
  UInt_t eventNumber, triggerTime, triggerTimeNs;
  Int_t firstAnt,secondAnt,maxAnt,corInd;
  Double_t deltaT,deltaTExpected;
  Double_t phiWave, phiMaxAnt;
  Double_t corPeak, corRMS;
  Double_t balloonLat, balloonLon, balloonAlt;
  Double_t heading,pitch,roll;
  

  Double_t thetaWave;


  int leftOpt, rightOpt;
  fGeomTool->getThetaPartners(par[0],leftOpt,rightOpt); 

  for(int i = 0; i <16; i++){
  deltaTArrayMod[i] = par[i];
  }

  //  deltaTArrayMod[rightOpt] = 0;

  for(entry=0;entry<numEntries;entry++) {
  

    corTree->GetEntry(entry);
    Long64_t headEntry=headTree->GetEntryNumberWithIndex(corSum->eventNumber);
    if(headEntry<0) 
      continue;
    headTree->GetEntry(headEntry);
     
    // if(header->triggerTimeNs*1e-9< 0.097 || header->triggerTimeNs*1e-9>0.1)
    if( (header->triggerTimeNs>0.3e6) || (header->triggerTimeNs<0.2e6) )  
      continue; 

    triggerTimeNs=header->triggerTimeNs;
    triggerTime=header->triggerTime;
    eventNumber=header->eventNumber;
  

    Long64_t bestEntry = adu5PatTree->GetEntryNumberWithBestIndex(header->triggerTime);
    if(bestEntry>-1) 
      adu5PatTree->GetEntry(bestEntry);
    else 
      continue;
      
    balloonLat=pat->latitude;
    balloonLon=pat->longitude;
    balloonAlt=pat->altitude;
    heading=pat->heading;
    pitch=pat->pitch;
    roll=pat->roll;
    UsefulAdu5Pat usefulPat(pat);
     
    for(corInd=0;corInd<19;corInd++) {


      //replace taylor dome
      //	 deltaTExpected=usefulPat.getDeltaTWillySeavey(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
      usefulPat.fSourceLongitude=0;
      deltaTExpected=usefulPat.getDeltaTTaylor(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);

      //	 histSimpleDtDiff->Fill(corSum->maxCorTimes[corInd]-deltaTExpected);
      int antInd=corSum->firstAnt[corInd]%16;
      //	 histDtDiffSep[antInd][labChip][corInd]->Fill(corSum->maxCorTimes[corInd]-deltaTExpected);

      firstAnt=corSum->firstAnt[corInd];
      secondAnt=corSum->secondAnt[corInd];
      deltaT=corSum->maxCorTimes[corInd];
      maxAnt=corSum->centreAntenna;
      phiWave=usefulPat.getPhiWave()*TMath::RadToDeg();
      phiMaxAnt=fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna)*TMath::RadToDeg();
      corPeak=corSum->maxCorVals[corInd];
      corRMS=corSum->rmsCorVals[corInd];
	
 if((deltaT - deltaTExpected)*(deltaT - deltaTExpected) < 1){
      phiAngle[0].push_back(phiWave);
      deltaTVec[0].push_back(deltaT - deltaTExpected + deltaTArrayMod[firstAnt] - deltaTArrayMod[secondAnt]);
      firstAntVec[0].push_back(firstAnt);
      secondAntVec[0].push_back(secondAnt);
 }
    }

    counter++; 
    
  }


  double deltaTArray[32][5000];
  double phiAngleArray[32][5000];

  double deltaTArrayCut[32][5000];
  double phiAngleArrayCut[32][5000];

  int middleAnt; 
  int leftAnt,rightAnt;

  int weirdOne[32];
  int meanPhi[32];
  int countArray[32];


  //fill arrays
  for(int ants = 0; ants < 32; ants++){
    int count = 0;
    int count2 = 0;
    double sumPhi = 0;
    bool true1 = false;
    bool true2 = false;

    fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
    
    for(int events = 1; events < phiAngle[0].size(); events++){
      int firstAntTemp = (int)firstAntVec[0][events];
      int secondAntTemp = (int)secondAntVec[0][events];
      int rightTemp = int(rightAnt);

      if((firstAntTemp == ants) &&  (secondAntTemp == rightTemp)){
	deltaTArray[ants][count] = deltaTVec[0][events];
	phiAngleArray[ants][count] = phiAngle[0][events];


	if(phiAngleArray[ants][count]<50){
	  true1 = true;
	}

	if(phiAngleArray[ants][count]>300){
	  true2=true;
	}

	count++;
      }

    }

    if(true1&&true2){
      weirdOne[ants]=1;
    }else{
      weirdOne[ants]=0;
    }

    countArray[ants] = count;

    for(int meanCalc = 0; meanCalc < count; meanCalc++){
      if(weirdOne[ants] && phiAngleArray[ants][meanCalc]> 200){
       	sumPhi = sumPhi + phiAngleArray[ants][meanCalc]-360;
      }else{
	sumPhi = sumPhi + phiAngleArray[ants][meanCalc];
      }

      count2++;
    }

    meanPhi[ants] = (sumPhi/count2);

    if(meanPhi[ants]<0){
      meanPhi[ants] = meanPhi[ants]+360;
    }

  }

  //make cuts
  for(int ants = 0; ants < 32; ants++){
    int count = 0;
    fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
    double sumPhi = 0;
 

    // cout << meanPhi[ants] << endl;

    for(int events = 0; events < countArray[ants]; events++){

      double lower  = meanPhi[ants]-20;
      double upper = meanPhi[ants]+20;

      if(lower > 360){
	lower = 0;
	upper = 20;
      }

      if(upper < 0){
	lower = 300;
	upper = 360;
      }

      if((phiAngleArray[ants][events] > lower ) && (phiAngleArray[ants][events]< upper)){
      
	phiAngleArrayCut[ants][count] = phiAngleArray[ants][events];
	deltaTArrayCut[ants][count] = deltaTArray[ants][events];
	count++;
      }

    }
   
    //cout << meanPhi[ants] << "  " <<  phiAngleArray[ants][10] << endl;

    TGraph *tempAntGraph  = new TGraph(count-1, phiAngleArrayCut[ants], deltaTArrayCut[ants]);
    
    //  canSurf->cd(ants+1);
    //canSurf->cd(1);
    tempAntGraph->SetMinimum(-1);
    tempAntGraph->SetMaximum(1);
    //tempAntGraph->Draw("ap");
    
    tempAntGraph->SetMarkerStyle(1);
   
    tempAntGraph->GetXaxis()->SetLimits(0,360);
if(ants <16){ 
    sumMean = sumMean + tempAntGraph->GetMean(2);
    sumMean2 = sumMean2 + tempAntGraph->GetMean(2)*tempAntGraph->GetMean(2);
 }
    //  cout << tempAntGraph->GetMean(2) << endl;


    theReturn = tempAntGraph->GetMean(2);

    delete tempAntGraph;    

  }

  delete event;
  delete hk; 
  delete header;
  delete pat;
  delete corSum;

  delete fpHead;
  delete fpGps;
  delete fpCor;

   cout << sumMean << "  " << sumMean2  << endl;

  return sumMean*sumMean + sumMean2;
  
}

