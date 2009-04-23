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
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "correlationTreeLoopOpt.C"

TCanvas *canSurf= new TCanvas("canSurf","canSurf");

//Double_t antOpt(Double_t *delta, TTree *headTree, TTree *adu5PatTree,TTree *corTree);
Double_t antOpt(Double_t *delta);

//Double_t antOpt(Double_t *delta, TTree *headTree, TTree *adu5PatTree,TTree *corTree){
Double_t antOpt(Double_t *delta){;

  TMultiGraph *myMG  = new TMultiGraph();

  Double_t histAdded = 0;
  Double_t histAdded2 = 0;

  //int middleAnt = delta[2];
  Double_t deltaR[16];  
  Double_t deltaZ[16];  
  Double_t deltaPhi[16];  
  Double_t deltaHeading[1];

  vector<Int_t> firstAntVec;
  vector<Int_t> secondAntVec;
  vector<double> deltaTVec;
  vector<Int_t> labChipVec;

  for(int r = 0; r <16; r++){
    deltaR[r] = delta[r];
  }
  
  for(int z = 0; z <16; z++){
    deltaZ[z] = 0;
  }
  
  for(int phi = 0; phi <16; phi++){
    //deltaPhi[phi] = 0;
    deltaPhi[phi] = delta[phi+16];
  }

  deltaHeading[0] = delta[18];

  // correlationTreeLoopOpt(18,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/","../../outFiles/","../../outFiles/",deltaR,deltaZ,deltaPhi,deltaHeading);

  //correlationTreeLoopOpt(18,"/unix/anita3/flight0809/root/","../outFiles/","../outFiles/",deltaR,deltaZ,deltaPhi,deltaHeading);
  correlationTreeLoopOpt(18,"/Users/simonbevan/Desktop/","../../outFiles/","../../outFiles/",deltaR,deltaZ,deltaPhi,deltaHeading);
  int run = 18;

  //char *baseDir = "http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/";
  char *baseDir = "/Users/simonbevan/Desktop/";
  // char *baseDir = "/unix/anita3/flight0809/root/";
  char *corTreeDir = "../../outFiles/";
  char *outputDir = "../../outFiles/";

  char eventName[100];
  char headerName[100];
  char hkName[100];
  char gpsName[100];
  char corrName[100];
  char outName[100];
  char fpTreeFile[100];

  //sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
  sprintf(eventName,"%s/run%d/eventFile%d.root",baseDir,run,run);
  sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
  sprintf(gpsName,"%s/run%d/gpsFile%d.root",baseDir,run,run);
  sprintf(corrName,"%s/corRun%d.root",corTreeDir,run);
  sprintf(outName,"%s/opt%d.root",outputDir,run);
  sprintf(fpTreeFile,"../../outFiles/deltaTFile18.root",run);

  TFile *fpTTree = new TFile(fpTreeFile);
  TTree *deltaTTree = (TTree*) fpTTree->Get("deltaTTree");


  //   TFile *fpHead = TFile::Open(headerName);
  //  TTree *headTree = (TTree*) fpHead->Get("headTree");
      
  //  TFile *fpGps = TFile::Open(gpsName);
  //  TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
   
  //  TFile *fpCor = TFile::Open(corrName);
  //  TTree *corTree = (TTree*) fpCor->Get("corTree");
  

  // RawAnitaHeader *header =0;
  // CorrelationSummary *corSum =0;
  // Adu5Pat *pat = 0;
  
  Int_t labChip;

  int leftAnt;
  int rightAnt;

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  //fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 
  double minValue = 0;

  // headTree->SetBranchAddress("header",&header);
  // headTree->BuildIndex("eventNumber");
  // adu5PatTree->BuildIndex("realTime");
  // adu5PatTree->SetBranchAddress("pat",&pat);
  // corTree->SetBranchAddress("cor",&corSum);
  // corTree->SetBranchAddress("labChip",&labChip);

  char plotCond[180];
  char plotTitle[180];
  char histNameOpt[180];

  TH1F *histAttenDiffOpt[16][1] = {{0}};
  for(int chip=0; chip<1;chip++){
    for(int atten = 0; atten <16; atten++){
      sprintf(plotTitle,"ANT %d - ANT %d",atten,atten+9);
      sprintf(histNameOpt,"histAttenDiffOpt%d_%d_%d",atten,atten+1,chip);
      histAttenDiffOpt[atten][chip]=new TH1F(histNameOpt,plotTitle,200,-1,1);
    }
  }

  for(int chip=0; chip<1;chip++){
    for(int ant=0;ant<16;ant++) {
      int middleAnt=ant; 
      int leftAnt,rightAnt;
      fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 

      sprintf(plotCond,"((firstAnt==%d && secondAnt==%d)) ",middleAnt,rightAnt,chip);	 
      sprintf(plotTitle,"Ant %d -  Ant %d",middleAnt,rightAnt);
      sprintf(histNameOpt,"histDtOpt_%d_%d",ant,chip);

      TH1F *histDt11Opt = new TH1F(histNameOpt,plotTitle,200,-1,1);
       
      deltaTTree->Project(histNameOpt,"deltaT-deltaTExpectedOpt",plotCond);
      histAttenDiffOpt[ant][chip]->Add(histDt11Opt);
      delete histDt11Opt;
    }
  }



  histAdded = 0;
  histAdded2 = 0;
  bool histOK = false;

  double sumMean = 0;

  for(int ants=0;ants<16;ants++) {
    for(int chip=0; chip<1;chip++){
      //histAdded = histAdded + histAttenDiffOpt[atten][chip]->GetMean()*histAttenDiffOpt[atten][chip]->GetMean();
      histAttenDiffOpt[ants][chip]->SetMinimum(-1.5);
      histAttenDiffOpt[ants][chip]->SetMaximum(1.5);
      histAttenDiffOpt[ants][chip]->Draw("ap");
    
      histAttenDiffOpt[ants][chip]->SetMarkerStyle(1);
   
      histAttenDiffOpt[ants][chip]->GetXaxis()->SetLimits(0,360);

      //cout << ants <<"  " <<histAttenDiffOpt[ant][chip]->GetMean(2) << endl;
      sumMean = sumMean + histAttenDiffOpt[ants][chip]->GetMean(2);

      if(ants == 8 || ants == 16 || ants == 12 || ants == 24){
	histAttenDiffOpt[ants][chip]->SetMarkerColor(8);
      }

      if(ants == 3 || ants == 7 || ants == 23 || ants == 31){
	histAttenDiffOpt[ants][chip]->SetMarkerColor(1);
      }
    
      if(ants == 9 || ants == 13 || ants == 18 || ants == 26){
	histAttenDiffOpt[ants][chip]->SetMarkerColor(2);
      }
    
    
      if(ants == 10 || ants == 14 || ants == 20 || ants == 28){
	histAttenDiffOpt[ants][chip]->SetMarkerColor(3);
      }
    
      if(ants == 11 || ants == 15 || ants == 22 || ants == 30){
	histAttenDiffOpt[ants][chip]->SetMarkerColor(4);
      }
    
      if(ants == 0 || ants == 4 || ants == 17 || ants == 25){
	histAttenDiffOpt[ants][chip]->SetMarkerColor(5);
      }


      if(ants == 2 || ants == 6 || ants == 21 || ants == 29){
	histAttenDiffOpt[ants][chip]->SetMarkerColor(6);
      }

      if(ants == 1 || ants == 5 || ants == 19 || ants == 27){
	histAttenDiffOpt[ants][chip]->SetMarkerColor(7);
      }

 

      

      //    if(ants==10){
      //      histAttenDiffOpt[ant][chip]->SetMarkerColor(3);
      //    }else{
      //      histAttenDiffOpt[ant][chip]->SetMarkerColor(ants);
      //    }
    
    
      //  myMG->Add(histAttenDiffOpt[ants][chip]);
      //myMG->Draw("p");


      histAdded = histAdded + histAttenDiffOpt[ants][chip]->GetMean();
      histAdded2 = histAdded + histAttenDiffOpt[ants][chip]->GetRMS();
      histOK = true;
    }
  }
 

  if(histOK){
    minValue = histAdded*histAdded+ histAdded2*histAdded2;
    //minValue = histAdded;
  }else{
    minValue = 10000;
    cout << "here" << endl;
  }

  //   cout << histAdded << endl;
  

  canSurf->Update();

  // delete header;
  // delete corSum;
  // delete pat;
  delete fGeomTool;

  //   delete deltaTTree;

  for(int i = 0; i<16; i++){
    delete histAttenDiffOpt[i][0];
  }


  delete fpTTree;
  //delete firstAntVec;
  //delete secondAntVec;
  //delete labChipVec;
  //delete deltaTVec;

  //    cout << "deltaR  " << deltaR[0] << "  "  << deltaHeading[0] << "  minValue  " << minValue << std::endl;
  cout << "The sum of the rings  " <<  histAdded  << " RMS  " << histAdded2<< " sumMean "<< sumMean <<endl;
  return minValue;
}

