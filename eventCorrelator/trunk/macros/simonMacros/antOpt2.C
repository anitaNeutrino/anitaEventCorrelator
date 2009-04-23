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

 Double_t deltaR[16];  
  Double_t deltaZ[16];  
  Double_t deltaPhi[16];  
  Double_t deltaHeading[1];

TCanvas *canSurf= new TCanvas("canSurf","canSurf");

void setAnt();
void setAnt(){
  for(int i = 0; i < 16; i++){
    deltaR[i] = 0;  
   deltaZ[i] = 0;  
    deltaPhi[i] = 0;  
    deltaHeading[i] = 0;
  }
}

void setAntValue(int x, double r, double phi, double z);
void setAntValue(int x, double r, double phi, double z){
  
    deltaR[x] = r;  
    deltaZ[x] = z;  
    deltaPhi[x] = phi;  
     deltaHeading[x] = 0;

  
}

void setArray();
void setArray(){
  for(int i = 0; i < 32; i++){

deltaTArrayMod[9] = -0.0222754;
deltaTArrayMod[1] = -0.0958421;
deltaTArrayMod[10] = 0.0649919;
deltaTArrayMod[2] = -0.183055;
deltaTArrayMod[11] = -0.176371;
deltaTArrayMod[3] = -0.0998352;
deltaTArrayMod[12] = 0.0986996;
deltaTArrayMod[4] = -0.133944;
deltaTArrayMod[13] = -0.164808;
deltaTArrayMod[5] = -0.029163;
deltaTArrayMod[14] = 0.0589837;
deltaTArrayMod[6] = -0.0364826;
deltaTArrayMod[15] = -0.0275802;
deltaTArrayMod[7] = -0.000186161;
deltaTArrayMod[8] = 0.24683;
deltaTArrayMod[0] = 0.000140377;

  }
}


void setValue(int x, double t);
void setValue(int x, double t){
  
  deltaTArrayMod[x] = t;
  
}

void printArray();
void printArray(){
  
  for(int i = 0; i < 16; i++){
    cout <<  deltaTArrayMod[i] << endl;;
  }
  cout << " " << endl;
}

void prinAnttArray();
void printAntArray(){
  
  for(int i = 0; i < 16; i++){
    cout <<  deltaR[i] << "  " << deltaZ[i] << "  "<< deltaPhi[i]<< endl;
  }
  cout << " " << endl;
}



double antOpt2(double *par);  
double antOpt2(double *par){
  


// for(int r = 0; r <16; r++){
//     deltaR[r] = par[r+5];
//   }
  
//   for(int z = 0; z <16; z++){
//     deltaZ[z] = 0;
//   }
  
  for(int phi = 0; phi <16; phi++){
    deltaPhi[phi] = par[phi+5];
  }

  deltaHeading[0] = 0;

  //correlationTreeLoopOpt(18,"/Users/simonbevan/Desktop/","../../outFiles/","../../outFiles/",deltaR,deltaZ,deltaPhi,deltaHeading);

  //  double deltaTArrayMod[32]={0};

  double theReturn = 0;
  double theReturn2 = 0;
  double sumMean = 0;
  double sumRMS = 0;

  TMultiGraph *myMG  = new TMultiGraph();

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corrName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  char baseDir[FILENAME_MAX];
  char *corTreeDir = "../../Outfiles";

  int run = 18;

  //sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
  sprintf(baseDir,"/Users/simonbevan/Desktop/");
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
  Double_t deltaT,deltaTExpected,deltaTExpectedOpt;
  Double_t phiWave, phiMaxAnt;
  Double_t corPeak, corRMS;
  Double_t balloonLat, balloonLon, balloonAlt;
  Double_t heading,pitch,roll;
  

  Double_t thetaWave;




  int leftOpt, rightOpt;
  fGeomTool->getThetaPartners(par[0],leftOpt,rightOpt); 

 //  deltaR[rightOpt] = par[1];  
//   deltaZ[rightOpt] = par[2];  
//   deltaPhi[rightOpt] = par[3];  
//   deltaTArrayMod[rightOpt] = par[4]; 
  

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
      deltaTExpectedOpt=usefulPat.getDeltaTTaylorOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd],deltaR,deltaZ,deltaPhi);

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
	deltaTVec[0].push_back(deltaT - deltaTExpectedOpt + deltaTArrayMod[firstAnt] - deltaTArrayMod[secondAnt]);
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
  //for(int ants = par[0]; ants < par[0]+1; ants++){
      for(int ants = 0; ants < 16; ants++){

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
      //       for(int ants = par[0]; ants < par[0]+1; ants++){
         for(int ants = 0; ants < 16; ants++){
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
    


    //  cout << tempAntGraph->GetMean(2) << endl;


    theReturn = theReturn + tempAntGraph->GetMean(2)*tempAntGraph->GetMean(2);
    theReturn2 = theReturn2 + tempAntGraph->GetRMS(2);

    // cout << theReturn << "  " << theReturn2 << endl;
    

    delete tempAntGraph;    

  }





  //make cuts
  for(int ants = 0; ants < 16; ants++){
    int count = 0;
    fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
    double sumPhi = 0;
 
    //cout <<  ants << "  " << meanPhi[ants] << endl;
    //  cout << [ants] << endl;

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

    TGraph *tempAntGraph  = new TGraph(count-1, phiAngleArrayCut[ants], deltaTArrayCut[ants]);
    
    if(ants <16){  
      //  canSurf->cd(ants+1);
      //canSurf->cd(1);
      tempAntGraph->SetMinimum(-1.5);
      tempAntGraph->SetMaximum(1.5);
      tempAntGraph->Draw("ap");
    
      tempAntGraph->SetMarkerStyle(1);
   
      tempAntGraph->GetXaxis()->SetLimits(0,360);
    
      //cout << ants <<"  " <<tempAntGraph->GetMean(2) << endl;
      sumMean = sumMean + tempAntGraph->GetMean(2);
      sumRMS = sumRMS + tempAntGraph->GetRMS(2);

      if(ants == 8 || ants == 16 || ants == 12 || ants == 24){
	tempAntGraph->SetMarkerColor(8);
      }

      if(ants == 3 || ants == 7 || ants == 23 || ants == 31){
	tempAntGraph->SetMarkerColor(1);
      }
    
      if(ants == 9 || ants == 13 || ants == 18 || ants == 26){
	tempAntGraph->SetMarkerColor(2);
      }
    
    
      if(ants == 10 || ants == 14 || ants == 20 || ants == 28){
	tempAntGraph->SetMarkerColor(3);
      }
    
      if(ants == 11 || ants == 15 || ants == 22 || ants == 30){
	tempAntGraph->SetMarkerColor(4);
      }
    
      if(ants == 0 || ants == 4 || ants == 17 || ants == 25){
	tempAntGraph->SetMarkerColor(5);
      }


      if(ants == 2 || ants == 6 || ants == 21 || ants == 29){
	tempAntGraph->SetMarkerColor(6);
      }

      if(ants == 1 || ants == 5 || ants == 19 || ants == 27){
	tempAntGraph->SetMarkerColor(7);
      }
    
      myMG->Add(tempAntGraph);
      myMG->Draw("p");

    }
  }


  canSurf->Update();

  cout <<"ant " <<  par[0] << " sum mean " << sumMean << " sumRMS " << sumRMS  << endl;

  

  delete event;
  delete hk; 
  delete header;
  delete pat;
  delete corSum;

  delete fpHead;
  delete fpGps ;
  delete fpCor;

  // for(int i = 0; i < 16; i++){
  // cout << theReturn*theReturn << "  " << deltaTArrayMod[i]  << endl;
  //}
  //cout << " " << endl;

  return sumMean*sumMean + sumRMS;
  
}

