#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "UsefulAdu5Pat.h"
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

void plotPrettyThings(int run, Long64_t entry, int ant);

void plotPrettyThingsEventNumber(int run, int eventNumber, int ant);


void plotPrettyThingsEventNumber(int run, int eventNumber, int ant)
{
  char headerName[FILENAME_MAX];    
  sprintf(headerName,"/unix/anita1/newRootData/run%d/timedHeadFile%d.root",run,run);
  TFile fp(headerName);
  TTree *headTree = (TTree*) fp.Get("headTree");
  headTree->BuildIndex("eventNumber");
  Long64_t entry=headTree->GetEntryNumberWithIndex(eventNumber);
  if(entry<0) {
    std::cout << "Couldn't find event " << eventNumber << " in run " << run << "\n";
    return;
  }
  plotPrettyThings(run,entry,ant);
}
  
void plotPrettyThings(int run, Long64_t entry, int ant) {

  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  sprintf(eventName,"/unix/anita1/newRootData/run%d/eventFile%d*.root",run,run);
  sprintf(headerName,"/unix/anita1/newRootData/run%d/timedHeadFile%d.root",run,run);
  sprintf(hkName,"/unix/anita1/newRootData/run%d/prettyHkFile%d.root",run,run);
  sprintf(gpsName,"/unix/anita1/newRootData/run%d/gpsFile%d.root",run,run);

  RawAnitaEvent *event = 0;
  TimedAnitaHeader *header =0;
  PrettyAnitaHk *hk = 0;
  Adu5Pat *pat =0;


  TChain *eventChain = new TChain("eventTree");
  eventChain->Add(eventName);
  eventChain->SetBranchAddress("event",&event);

  TFile *fpHead = new TFile(headerName);
  TTree *headTree = (TTree*) fpHead->Get("headTree");
  headTree->SetBranchAddress("header",&header);

  TFile *fpHk = new TFile(hkName);
  TTree *prettyHkTree = (TTree*) fpHk->Get("prettyHkTree");
  prettyHkTree->SetBranchAddress("hk",&hk);

  TFile *fpGps = new TFile(gpsName);
  TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
  adu5PatTree->SetBranchAddress("pat",&pat);


  //Friends only seem to work with TTree::Draw and similar commands
  //if you are manually calling GetEntry (i.e in a loop) you must call
  //the GetEntry for each tree separately.
  //  eventChain->AddFriend(headTree);
  //  eventChain->AddFriend(prettyHkTree);
    
  //Stupidly most do this to be perfectly safe  
  eventChain->GetEntry(entry);
  headTree->GetEntry(entry);
  prettyHkTree->GetEntry(entry);
  adu5PatTree->GetEntry(entry);
    
  
  PrettyAnitaEvent realEvent(event,WaveCalType::kVTFullJWPlusFastClockZero,hk);
  UsefulAdu5Pat usefulPat(pat);
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
  TCanvas *canWave = realEvent.getTenWaveformCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canPower = realEvent.getSixPowerEnvelopeCanvas(ant, AnitaPol::kVertical);
  TCanvas *canInter = realEvent.getSixInterpolatedCanvas(ant, AnitaPol::kVertical);
  canInter->SetWindowSize(1200,400);


  //  TCanvas *canFFT = realEvent.getSixFFTPowerCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canCor = realEvent.getSixCorrelationCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canElevenCor = realEvent.getElevenCorrelationCanvas(ant, AnitaPol::kVertical);
  TCanvas *canElevenIntCor = realEvent.getElevenInterpolationCorrelationCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canIntCor = realEvent.getSixInterpolatedCorrelationCanvas(ant, AnitaPol::kVertical);
  // 
//   CorrelationSummary *corSum=realEvent.getCorrelationSummary(-1,AnitaPol::kVertical,1./(2.6*16.));
//   //CorrelationSummary *corSum=realEvent.getCorrelationSummary(-1,AnitaPol::kVertical             );
//   cout << "\n\n\n";
//   for(int corInd=0;corInd<19;corInd++) {
//      Double_t deltaTExpected=usefulPat.getDeltaTWillySeavey(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
//      cout << corSum->firstAnt[corInd] << "\t" << corSum->secondAnt[corInd] << "\t" << deltaTExpected << "\t" << corSum->maxCorTimes[corInd] << "\t" << (corSum->maxCorTimes[corInd]-deltaTExpected) << "\n";
//   }
}

