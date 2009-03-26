#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "RawAnitaHeader.h"
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
#include "TStopwatch.h"
#include <iostream>
#include <fstream>

void makeCorrelationRunTree(int run, int numEnts, char *baseDir, char *outDir=0);

int main(int argc, char **argv) {
  int run=17;
  if(argc>1) {
    run=atoi(argv[1]);
  }
  int numEnts=0;
  if(argc>2) {
    numEnts=atoi(argv[2]);
  }
  std::cout << "Making correlation summary tree for run: " << run << "\n";

  // makeCorrelationRunTree(run,0,"/Users/simonbevan/Desktop/","/Users/simonbevan/ANITA/outfiles/");

  //makeCorrelationRunTree(run,0,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/","/home/rjn/anita/data/corTrees");
  // makeCorrelationRunTree(run,0,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/","../outFiles/");

  TStopwatch stopy;
  stopy.Start();
  makeCorrelationRunTree(run,0,"/Users/simonbevan/Desktop/","/Users/simonbevan/ANITA/outfiles/");
  //  makeCorrelationRunTree(run,numEnts,"/unix/anita3/flight0809/root","/unix/anita3/rjn/corTreesSlow");
  stopy.Stop();
  std::cout << "Run " << run << "\t" << numEnts << " events \n";
  std::cout << "CPU Time: " << stopy.CpuTime() << "\t" << "Real Time: "
	    << stopy.RealTime() << "\n";


  
   //makeCorrelationRunTree(run,0,"/unix/anita3/flight0809/root","/unix/anita3/rjn/corTrees");

}
  
void makeCorrelationRunTree(int run, int numEnts, char *baseDir, char *outDir) {
  
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  
  //The locations of the event, header and gps files 
  // The * in the evnt file name is a wildcard for any _X files  
  sprintf(eventName,"%s/run%d/calEventFile%d*.root",baseDir,run,run);
  sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
  sprintf(gpsName,"%s/run%d/gpsEvent%d.root",baseDir,run,run);

  Int_t useCalibratedFiles=0;

  //Define and zero the class pointers
  RawAnitaEvent *event = 0;
  CalibratedAnitaEvent *calEvent = 0;
  RawAnitaHeader *header =0;
  Adu5Pat *pat = 0;
  
  //Need a TChain for the event files as they are split in to sub files to ensure no single file is over 2GB
  TChain *eventChain = new TChain("eventTree");
  eventChain->Add(eventName);

  //Here we check if there are any entries in the tree of CalibratedAnitaEvent objects
  if(eventChain->GetEntries()>0) {
    //If there are entries we can set the branc address
    eventChain->SetBranchAddress("event",&calEvent);
    useCalibratedFiles=1;
  }
  else {
    //If there aren't any entries we can try raw event files instead
    sprintf(eventName,"%s/run%d/eventFile%d*.root",baseDir,run,run);
    eventChain->Add(eventName);
    eventChain->SetBranchAddress("event",&event);
  }
  
  //Now open the header and GPS files
  TFile *fpHead = TFile::Open(headerName);
  TTree *headTree = (TTree*) fpHead->Get("headTree");
  headTree->SetBranchAddress("header",&header);

  TFile *fpGps = TFile::Open(gpsName);
  TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
  adu5PatTree->SetBranchAddress("pat",&pat);
  
  //The index is necessary until we have an interpolated file
  //adu5PatTree->BuildIndex("realTime");

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
  Double_t thetaWave,phiWave;
  Int_t labChip;
  TTree *corTree = new TTree("corTree","Tree of Correlation Summaries");
  corTree->Branch("cor","CorrelationSummary",&theCor);
  corTree->Branch("thetaWave",&thetaWave,"thetaWave/D");
  corTree->Branch("phiWave",&phiWave,"phiWave/D");
  corTree->Branch("labChip",&labChip,"labChip/I");


  std::cout << headTree->GetEntries() << "  " << eventChain->GetEntries() << std::endl;


  Long64_t maxEntry=headTree->GetEntries(); 
  if(numEnts && maxEntry>numEnts) maxEntry=numEnts;

  Int_t starEvery=maxEntry/10;
  if(starEvery==0) starEvery=1;
  
  std::cout <<  "There are " << maxEntry << " events to proces\n";
  for(Long64_t entry=0;entry<maxEntry;entry++) {
     if(entry%starEvery==0) std::cerr << "*";

     //Get header
     headTree->GetEntry(entry);

     //    if( (header->triggerTimeNs>0.4e6) || (header->triggerTimeNs<0.25e6) )  
     //Now cut to only process the Taylor Dome pulses
          if( (header->triggerTimeNs>0.5e6) || (header->triggerTimeNs<0.2e6) )  
       continue; 
     
     //Get event
     eventChain->GetEntry(entry);

     //Get GPS by finding the entry closest to the trigger time
     //Long64_t bestEntry = adu5PatTree->GetEntryNumberWithBestIndex(header->triggerTime);
     //if(bestEntry>-1) 
       adu5PatTree->GetEntry(entry);
       //else {
       //std::cerr << "No GPS for event " << header->eventNumber << "\n";
       //continue;
       // }

     //Now we can make a PrettyAnitaEvent (or a UsefulAnitaEvent)
     PrettyAnitaEvent *realEvent=0;
     if(useCalibratedFiles) {
       //If we have CalibratedAnitaEvent then the constructor is just
       realEvent = new PrettyAnitaEvent(calEvent);
     }
     else {
       //If we have RawAnitaEvent then we have to specify the calibration option
       realEvent = new PrettyAnitaEvent(event,WaveCalType::kVTFullAGCrossCorClock,header);
     }
     labChip=realEvent->getLabChip(1);
     
     // UsefulAdu5Pat contains some generically useful GPS orientation thingies
     UsefulAdu5Pat usefulPat(pat);
     usefulPat.getThetaAndPhiWaveWillySeavey(thetaWave,phiWave);
     int ant=realEvent->getMaxAntenna(AnitaPol::kVertical);
     //     int ant=fGeomTool->getUpperAntNearestPhiWave(phiWave);
     Double_t deltaT= 1. / (2.6*16);
     //     std::cout << deltaT << std::endl;
     //     if(entry%100==0)
     //       std::cout << ant << "\t" << realEvent.getMaxAntenna(AnitaPol::kVertical) << std::endl;
     theCor =realEvent->getCorrelationSummary(ant,AnitaPol::kVertical,deltaT);
     corTree->Fill();     
     delete theCor;
     if(realEvent) delete realEvent;
  }
  std::cerr << "\n";
  corTree->AutoSave();
  fpOut->Close();
}
     
