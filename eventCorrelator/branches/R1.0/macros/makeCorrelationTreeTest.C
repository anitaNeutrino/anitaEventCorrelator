#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "CorrelationSummary.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>

void makeCorrelationTreeTest(int run, int entry, int ant);

  
void makeCorrelationTreeTest(int run, int entry, int ant) {

  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  sprintf(eventName,"/unix/anita1/webData/firstDay/run%d/eventFile%d*.root",run,run);
  sprintf(headerName,"/unix/anita1/webData/firstDay/run%d/timedHeadFile%d.root",run,run);
  sprintf(hkName,"/unix/anita1/webData/firstDay/run%d/prettyHkFile%d.root",run,run);

  RawAnitaEvent *event = 0;
  TimedAnitaHeader *header =0;
  PrettyAnitaHk *hk = 0;
  
  TChain *eventChain = new TChain("eventTree");
  eventChain->Add(eventName);
  eventChain->SetBranchAddress("event",&event);

  TFile *fpHead = new TFile(headerName);
  TTree *headTree = (TTree*) fpHead->Get("headTree");
  headTree->SetBranchAddress("header",&header);

  TFile *fpHk = new TFile(hkName);
  TTree *prettyHkTree = (TTree*) fpHk->Get("prettyHkTree");
  prettyHkTree->SetBranchAddress("hk",&hk);


  //Friends only seem to work with TTree::Draw and similar commands
  //if you are manually calling GetEntry (i.e in a loop) you must call
  //the GetEntry for each tree separately.
  //  eventChain->AddFriend(headTree);
  //  eventChain->AddFriend(prettyHkTree);
    
  //Stupidly most do this to be perfectly safe  
  eventChain->GetEntry(entry);
  headTree->GetEntry(entry);
  prettyHkTree->GetEntry(entry);
    
  
  PrettyAnitaEvent realEvent(event,WaveCalType::kVTFullJWPlus,hk);
  cout << realEvent.eventNumber << " " << header->eventNumber << endl;
//   TGraph *gr = realEvent.getGraph(14,AnitaPol::kVertical);
//   Double_t x,y;
//   Double_t lastx=0;
  
//   for(int i=0;i<gr->GetN();i++) {
//     gr->GetPoint(i,x,y);
//     cout << i << "\t" << x-lastx << endl;
//     lastx=x;
//   }
 

//  TCanvas *canWave = realEvent.getSixWaveformCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canPower = realEvent.getSixPowerEnvelopeCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canInter = realEvent.getSixInterpolatedCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canFFT = realEvent.getSixFFTPowerCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canCor = realEvent.getSixCorrelationCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canElevenCor = realEvent.getElevenCorrelationCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canElevenIntCor = realEvent.getElevenInterpolationCorrelationCanvas(ant, AnitaPol::kVertical);
 // TCanvas *canIntCor = realEvent.getSixInterpolatedCorrelationCanvas(ant, AnitaPol::kVertical);
  // 
  CorrelationSummary *theCor=0;
  TFile *fp = new TFile("corTest.root","RECREATE");
  TTree *corTree = new TTree("corTree","Tree of Correlation Summaries");
  corTree->Branch("cor","CorrelationSummary",&theCor);

  
  //  theCor=realEvent.getCorrelationSummary(AnitaPol::kVertical);
  //  corTree->Fill();
  //  delete theCor;

  //  for(Double_t deltaT=1./2.6; deltaT>1./(2.6*32); deltaT/=2.6) {
  for(double i=1;i<32;i+=1) {
     Double_t deltaT= 1. / (2.6*i);
     //     cout << deltaT << endl;
     theCor =realEvent.getCorrelationSummary(AnitaPol::kVertical,deltaT);
     corTree->Fill();     
     delete theCor;
  }
  corTree->AutoSave();
  fp->Close();
}

