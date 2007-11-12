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

  
void makeCorrelationRunTree(int run, int numEnts, char *outDir) {
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();

  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  sprintf(eventName,"/unix/anita1/webData/firstDay/run%d/eventFile%d*.root",run,run);
  sprintf(headerName,"/unix/anita1/webData/firstDay/run%d/timedHeadFile%d.root",run,run);
  sprintf(hkName,"/unix/anita1/webData/firstDay/run%d/prettyHkFile%d.root",run,run);
  sprintf(gpsName,"/unix/anita1/webData/firstDay/run%d/gpsFile%d.root",run,run);

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
  cout << outName << endl;
  fpOut= new TFile(outName,"RECREATE");
  TTree *corTree = new TTree("corTree","Tree of Correlation Summaries");
  corTree->Branch("cor","CorrelationSummary",&theCor);

  Long64_t maxEntry=headTree->GetEntries(); 
  if(numEnts && maxEntry>numEnts) maxEntry=numEnts;

  Int_t starEvery=maxEntry/1000;
  if(starEvery==0) starEvery=1;
  
  Double_t thetaWave,phiWave;

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
     
     
     PrettyAnitaEvent realEvent(event,WaveCalType::kVTFullJWPlusClockZero,hk);
     UsefulAdu5Pat usefulPat(pat);
     usefulPat.getThetaAndPhiWaveWillySeavey(thetaWave,phiWave);
     int ant=fGeomTool->getUpperAntNearestPhiWave(phiWave);
	 //     int ant=9;
     //     cout << realEvent.eventNumber << " " << header->eventNumber << endl;
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
     
     
     //  theCor=realEvent.getCorrelationSummary(AnitaPol::kVertical);
     //  corTree->Fill();
     //  delete theCor;
     
     //  for(Double_t deltaT=1./2.6; deltaT>1./(2.6*32); deltaT/=2.6) {
     //     for(double i=1;i<32;i+=0.1) {
     Double_t deltaT= 1. / (2.6*24.);
     //     cout << deltaT << endl;
     //     if(entry%100==0)
     //	std::cout << ant << "\t" << realEvent.getMaxAntenna(AnitaPol::kVertical) << std::endl;
     theCor =realEvent.getCorrelationSummary(ant,AnitaPol::kVertical,deltaT);
     corTree->Fill();     
     delete theCor;
  //    delete event;
//      delete pat;
//      delete header;
//      delete hk;
     
  }
  std::cerr << "\n";
  corTree->AutoSave();
  fpOut->Close();
}
     
