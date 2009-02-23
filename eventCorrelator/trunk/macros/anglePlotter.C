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



void anglePlotter();

  
void anglePlotter(){
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


  double deltaTArrayMod[32];
//   deltaTArrayMod[9] =  -0.00415702;
//   deltaTArrayMod[1] =  0.0442874;
//   deltaTArrayMod[10] =  -0.249513;
//   deltaTArrayMod[2] =  0.493758;
//   deltaTArrayMod[11] =  -0.550521;
//   deltaTArrayMod[3] =  0.501072;
//   deltaTArrayMod[12] =  -0.732701;
//   deltaTArrayMod[4] =  0.993376;
//   deltaTArrayMod[13] =  -0.95656;
//   deltaTArrayMod[5] =  0.862654;
//   deltaTArrayMod[14] =  -0.911761;
//   deltaTArrayMod[6] =  1.02809;
//   deltaTArrayMod[15] =  -0.982902;
//   deltaTArrayMod[7] =  0.932202;
//   deltaTArrayMod[8] =  -1.12693;


     deltaTArrayMod[0] = -0.0183779;
     deltaTArrayMod[1] = -0.144295;
     deltaTArrayMod[2] =0.139706;
     deltaTArrayMod[3] = -0.0228605;
     deltaTArrayMod[4] =0.298226;
     deltaTArrayMod[5] = 0.00239688;
     deltaTArrayMod[6] =-0.000770198;
     deltaTArrayMod[7] = -0.254639;
     deltaTArrayMod[8] =0.143753;
     deltaTArrayMod[9] = 0.10019;
     deltaTArrayMod[10] =0.020126;
     deltaTArrayMod[11] = -0.11267;
     deltaTArrayMod[12] =-0.123496;
     deltaTArrayMod[13] = -0.179272;
     deltaTArrayMod[14] = 0.0338776;
     deltaTArrayMod[15]  = 0.117575;


  for(int i = 16; i < 32; i++){
    deltaTArrayMod[i] = 0;
  }


  double sumMean = 0;
  double sumMean2 = 0;

  // sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
  sprintf(baseDir,"/unix/anita3/flight0809/root/");
  sprintf(eventName,"%s/run%d/eventFile%d.root",baseDir,run,run);
  sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
  sprintf(gpsName,"%s/run%d/gpsFile%d.root",baseDir,run,run);
  sprintf(corrName,"%s/corRun%d.root",corTreeDir,run);

  TCanvas *canSurf= new TCanvas("canSurf","canSurf");
  canSurf->Divide(1,2);

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


  double deltaR[32]={0};
  deltaR[0] = -1.52347e-05;  
  // deltaR[1] =0.0469885;
   deltaR[1] =-0.039885;
  deltaR[2] =4.76685e-05;  
  deltaR[3] =-0.033546;    
  //deltaR[4] =0.258364;
  deltaR[5] =-0.000207312;    
  deltaR[6] =-0.0175568;   
  deltaR[7] =-2.78447e-06; 
  deltaR[8] =0.018968;    
  deltaR[9] =-4.6708e-07;  
  deltaR[10] =0.00530499;  
  deltaR[11] =5.82057e-05; 
  deltaR[12] =0.00285399 ; 
  deltaR[13] =0.00365397;  
  deltaR[14] =1.11022e-06; 
  deltaR[15] =0.00454388;  

  double deltaPhi[32]={0};
 deltaPhi[0] =  -0.000654457;
 deltaPhi[1] =  -0.00513332;
 deltaPhi[2] =  -5.3362e-06;
 deltaPhi[3] = -4.5843e-07;
 deltaPhi[4] =  -0.0073673;
 deltaPhi[5] = -6.96774e-07;
 deltaPhi[6] = -0.000760718;
 deltaPhi[7] =  0.00209363;
 deltaPhi[8] = -0.000678059;

double deltaZ[32]={0};

  //sprintf(eventName,"/unix/anita1/webData/firstDay/run%d/eventFile%d*.root",run,run);
  //sprintf(headerName,"/unix/anita1/webData/firstDay/run%d/timedHeadFile%d.root",run,run);
  //sprintf(hkName,"/unix/anita1/webData/firstDay/run%d/prettyHkFile%d.root",run,run);
  //sprintf(gpsName,"/unix/anita1/webData/firstDay/run%d/gpsFile%d.root",run,run);
  //sprintf(corrName,"/unix/anita1/rjn/corTree16/corRun%d.root",run);
  //sprintf(outName,"deltaTFileSlowClock%d.root",run);
   
  RawAnitaEvent *event = 0;
  PrettyAnitaHk *hk = 0;
   
  RawAnitaHeader *header =0;
  Adu5Pat *pat =0;
  CorrelationSummary *corSum =0;
   
  //  TChain *eventChain = new TChain("eventTree");
  // eventChain->Add(eventName);
  // eventChain->SetBranchAddress("event",&event);
   
  //  TFile *fpHk = new TFile(hkName);
  //    TTree *prettyHkTree = (TTree*) fpHk->Get("prettyHkTree");
  //    prettyHkTree->SetBranchAddress("hk",&hk);
   
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


 

  // int leftOpt, rightOpt;
  //  fGeomTool->getThetaPartners(h,leftOpt,rightOpt); 


  //   deltaTArrayMod[rightOpt] = 0;

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
      deltaTExpected=usefulPat.getDeltaTTaylorOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd],deltaR,deltaZ,deltaPhi);




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
	deltaTVec[0].push_back(deltaT - deltaTExpected + deltaTArrayMod[firstAnt] + deltaTArrayMod[secondAnt]);
	phiAngle[0].push_back(phiWave);      
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

  TMultiGraph *myMG  = new TMultiGraph();
  TMultiGraph *myMG2  = new TMultiGraph();


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

  //  //calculate mean
  //   for(int ants = 0; ants < 32; ants++){
  //     int count = 0;
  //     fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
  //     double sumPhi = 0;
  //     for(int events = 1; events < phiAngle[0].size(); events++){

  //       if(weirdOne[ants] && phiAngleArray[ants][events] < 100){
  //        	sumPhi = sumPhi + phiAngleArray[ants][events]+360;
  //       }else{
  // 	sumPhi = sumPhi + phiAngleArray[ants][events];
  //       }

  //       count++;

  //     }
  //     meanPhi[ants] = (sumPhi/count-1);

  

  //   }
  

  //make cuts
  for(int ants = 0; ants < 16; ants++){
    int count = 0;
    fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
    double sumPhi = 0;
 
cout <<  ants << "  " << meanPhi[ants] << endl;
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
   
    //cout << meanPhi[ants] << "  " <<  phiAngleArray[ants][10] << endl;

    TGraph *tempAntGraph  = new TGraph(count-1, phiAngleArrayCut[ants], deltaTArrayCut[ants]);
    
    if(ants <16){  
      //  canSurf->cd(ants+1);
      canSurf->cd(1);
      tempAntGraph->SetMinimum(-1.5);
      tempAntGraph->SetMaximum(1.5);
      tempAntGraph->Draw("ap");
    
      tempAntGraph->SetMarkerStyle(1);
   
      tempAntGraph->GetXaxis()->SetLimits(0,360);
    


      //cout << ants <<"  " <<tempAntGraph->GetMean(2) << endl;
      sumMean = sumMean + tempAntGraph->GetMean(2);

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

 



      //    if(ants==10){
      //      tempAntGraph->SetMarkerColor(3);
      //    }else{
      //      tempAntGraph->SetMarkerColor(ants);
      //    }
    
    
      myMG->Add(tempAntGraph);
      myMG->Draw("p");




    }else{

      canSurf->cd(2);
      
      tempAntGraph->SetMinimum(-1);
      tempAntGraph->SetMaximum(1);
      tempAntGraph->Draw("ap");
    
      tempAntGraph->SetMarkerStyle(1);
    
      //cout << tempAntGraph->GetMean(2) << endl;
      sumMean2 = sumMean2 + tempAntGraph->GetMean(2);

      if(ants == 8 || ants == 16 || ants == 12 || ants == 24){
	tempAntGraph->SetMarkerColor(0);
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





      //     if(ants-16 == 10){
      //       tempAntGraph->SetMarkerColor(3);
      //     }else{
      //       tempAntGraph->SetMarkerColor(ants-16);
      //     }

      tempAntGraph->GetXaxis()->SetLimits(0,360);
    
      // canSurf->Update();
    
      tempAntGraph->SetTitle("Surf Correlations");
      tempAntGraph->GetXaxis()->SetTitle("Phi (degrees)");
      tempAntGraph->GetYaxis()->SetTitle("actual - expected (ns)");

      myMG2->Add(tempAntGraph);
      myMG2->Draw("p");

      myMG2->Add(tempAntGraph);
    }
  }

  

  cout << "******" << endl;
  cout << sumMean << endl;
  cout << sumMean2 << endl;
  cout << sumMean2+sumMean << endl;
  cout << "******" << endl;

  //canSurf->Update();

  //  canSurf->Update();


}

