#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
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

void plotPrettyThings(int run, int entry, int ant);

  
void plotPrettyThings(int run, int entry, int ant) {

  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  sprintf(eventName,"/home/rjn/anita/data/run%d/eventFile%d*.root",run,run);
  sprintf(headerName,"/home/rjn/anita/data/run%d/timedHeadFile%d.root",run,run);
  sprintf(hkName,"/home/rjn/anita/data/run%d/prettyHkFile%d.root",run,run);

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
    
    
  PrettyAnitaEvent realEvent(event,WaveCalType::kVTFullJWPlusFastClockZero,hk);
  cout << realEvent.eventNumber << " " << header->eventNumber << endl;
  TCanvas *canWave = realEvent.getSixWaveformCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canPower = realEvent.getSixPowerEnvelopeCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canInter = realEvent.getSixInterpolatedCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canFFT = realEvent.getSixFFTPowerCanvas(ant, AnitaPol::kVertical);
 TCanvas *canCor = realEvent.getSixCorrelationCanvas(ant, AnitaPol::kVertical);
 TCanvas *canElevenCor = realEvent.getElevenCorrelationCanvas(ant, AnitaPol::kVertical);
 // TCanvas *canElevenIntCor = realEvent.getElevenInterpolationCorrelationCanvas(ant, AnitaPol::kVertical);
 // TCanvas *canIntCor = realEvent.getSixInterpolatedCorrelationCanvas(ant, AnitaPol::kVertical);
}

