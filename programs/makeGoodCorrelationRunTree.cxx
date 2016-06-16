#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "AnitaGeomTool.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
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

void makeCorrelationRunTree(int run, int numEnts=0, char *outDir=0);

int main(int argc, char **argv) {
  int run=1028;
  if(argc>1) {
    run=atoi(argv[1]);
  }
  std::cout << "Making correlation summary tree for run: " << run << "\n";
  makeCorrelationRunTree(run,0,"/unix/anita1/rjn/peakCut5");
}
  
void makeCorrelationRunTree(int run, int numEnts, char *outDir) {
  //AnitaGeomTool *fGeomTool = 
  //  AnitaGeomTool::Instance();

  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  sprintf(eventName,"/unix/anita1/newRootData/run%d/eventFile%d*.root",run,run);
  sprintf(headerName,"/unix/anita1/newRootData/run%d/timedHeadFile%d.root",run,run);
  sprintf(hkName,"/unix/anita1/newRootData/run%d/prettyHkFile%d.root",run,run);
  sprintf(gpsName,"/unix/anita1/newRootData/run%d/gpsFile%d.root",run,run);

  RawAnitaEvent *event = 0;
  TimedAnitaHeader *header =0;
  PrettyAnitaHk *hk = 0;
  Adu5Pat *pat = 0;
  
  TChain *eventChain = new TChain("eventTree");
  eventChain->Add(eventName);
  eventChain->SetBranchAddress("event",&event);

  TFile *fpHead = new TFile(headerName);
  TTree *headTree = (TTree*) fpHead->Get("headTree");
  headTree->SetBranchAddress("header",&header);

  TFile *fpGps = new TFile(gpsName);
  TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
  adu5PatTree->SetBranchAddress("pat",&pat);

  TFile *fpHk = new TFile(hkName);
  TTree *prettyHkTree = (TTree*) fpHk->Get("prettyHkTree");
  prettyHkTree->SetBranchAddress("hk",&hk);


  //Make output files
  CorrelationSummary *theCor=0;
  TFile *fpOut;
  if(outDir) {
     sprintf(outName,"%s/corRun%d.root",outDir,run);
  }
  else {
     sprintf(outName,"corRun%d.root",run);
  }
  std::cout << outName << std::endl;
  fpOut= new TFile(outName,"RECREATE");
  TTree *corTree = new TTree("corTree","Tree of Correlation Summaries");
  corTree->Branch("cor","CorrelationSummary",&theCor);

  Long64_t maxEntry=headTree->GetEntries(); 
  if(numEnts && maxEntry>numEnts) maxEntry=numEnts;

  Int_t starEvery=maxEntry/1000;
  if(starEvery==0) starEvery=1;
  
  Double_t thetaWave,phiWave;
  std::cout <<  "There are " << maxEntry << " events to process\n";
  for(Long64_t entry=0;entry<maxEntry;entry++) {
     if(entry%starEvery==0) std::cerr << "*";
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
     //     usefulPat.getThetaAndPhiWaveWillySeavey(thetaWave,phiWave);
     Double_t peak=0;
     int ant=realEvent.getMaxAntenna(AnitaPol::kVertical,&peak);
     //     int ant=fGeomTool->getUpperAntNearestPhiWave(phiWave);
     Double_t deltaT= 1. / (2.6*16);
     //     std::cout << deltaT << std::endl;
     //     if(entry%100==0)
     //     std::cout << entry << "\t" << ant << "\t" << peak << std::endl; 
     if(peak<5) continue;
     //       std::cout << ant << "\t" << realEvent.getMaxAntenna(AnitaPol::kVertical) << std::endl;
     theCor =realEvent.getCorrelationSummary(ant,AnitaPol::kVertical,deltaT);
     corTree->Fill();     
     delete theCor;
  }
  std::cerr << "\n";
  corTree->AutoSave();
  fpOut->Close();
}
     
