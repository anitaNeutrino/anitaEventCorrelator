#include "AnitaConventions.h"
#include "UsefulAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "AnitaGeomTool.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
#include "FFTtools.h"
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

void makeSingleEventCorrelation(int run, int seedEvent, const char *baseDir, const char *outDir=0);

int main(int argc, char **argv) {
  if(argc<2) {
    std::cout << "Usage " << argv[0] << " <run> <seedEvent> <baseDir> <outDir>\n";
    return -1;
  }
  int run=17;
  if(argc>1) {
    run=atoi(argv[1]);
  }
  int seedEvent=0;
  if(argc>2) {
    seedEvent=atoi(argv[2]);
  }
  std::cout << "Making single event correlation tree for run: " << run << "\t" << seedEvent << "\n";

  // makeSingleEventCorrelation(run,0,"/Users/simonbevan/Desktop/","/Users/simonbevan/ANITA/outfiles/");

  //makeSingleEventCorrelation(run,0,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/","/home/rjn/anita/data/corTrees");
  // makeSingleEventCorrelation(run,0,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/","../outFiles/");

  TStopwatch stopy;
  stopy.Start();
  //  makeSingleEventCorrelation(run,0,"/Users/simonbevan/Desktop/","/Users/simonbevan/ANITA/outfiles/");
  makeSingleEventCorrelation(run,seedEvent,"/Users/rjn/Dropbox/Projects/anita/anita3/data/root","/Users/rjn/Dropbox/Projects/anita/anita3/singleEvent");
  stopy.Stop();
  std::cout << "Run " << run << "\n";
  std::cout << "CPU Time: " << stopy.CpuTime() << "\t" << "Real Time: "
	    << stopy.RealTime() << "\n";


  
   //makeSingleEventCorrelation(run,0,"/unix/anita3/flight0809/root","/unix/anita3/rjn/corTrees");

}
  
void makeSingleEventCorrelation(int run, int seedEvent, const char *baseDir, const char *outDir) {
  
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  
  //The locations of the event, header and gps files 
  // The * in the evnt file name is a wildcard for any _X files  
  sprintf(eventName,"%s/run%d/calEventFile%d*.root",baseDir,run,run);
  sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);

  Int_t useCalibratedFiles=0;

  //Define and zero the class pointers
  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  RawAnitaEvent *event = 0;
  CalibratedAnitaEvent *calEvent = 0;
  RawAnitaHeader *header =0;
  
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


  headTree->BuildIndex("eventNumber");
  Long64_t seedEntry= headTree->GetEntryNumberWithIndex(seedEvent);
  std::cout << "\nSeed Entry: " << seedEntry << "\n";


  //Get header
  headTree->GetEntry(seedEntry);
  //Get event
  eventChain->GetEntry(seedEntry);
  
  
  //Make output files
  TFile *fpOut;
  if(outDir) {
    sprintf(outName,"%s/corEvent%d_%d.root",outDir,run,seedEvent);
  }
  else {
    sprintf(outName,"corEvent%d_%d.root",run,seedEvent);
  }
  std::cout << outName << std::endl;
  fpOut= new TFile(outName,"RECREATE");

  Double_t peakCorr[3][3];
  Double_t peakTime[3][3];
  Double_t peakToPeakBottomTop[3];
  Double_t peakToPeak[3][3];
  Double_t rms[3][3];
  Double_t sumPeakCorr;
  Double_t sumPeakCorrSq;
  Int_t eventNumber=0;
  
  
  TTree *corTree = new TTree("corTree","corTree");
  corTree->Branch("eventNumber",&eventNumber,"eventNumber/I");
  corTree->Branch("peakCorr",peakCorr,"peakCorr[3][3]/D");
  corTree->Branch("peakTime",peakTime,"peakTime[3][3]/D");
  corTree->Branch("peakToPeakBottomTop",peakToPeakBottomTop,"peakToPeakBottomTop[3]/D");
  corTree->Branch("peakToPeak",peakToPeak,"peakToPeak[3][3]/D");
  corTree->Branch("rms",rms,"rms[3][3]/D");  
  corTree->Branch("sumPeakCorr",&sumPeakCorr,"sumPeakCorr/D");
  corTree->Branch("sumPeakCorrSq",&sumPeakCorrSq,"sumPeakCorrSq/D");

  
  std::cout << headTree->GetEntries() << "  " << eventChain->GetEntries() << std::endl;

  //Now we can make a UsefulAnitaEvent (or a UsefulAnitaEvent)
  UsefulAnitaEvent *testEvent=0;
  if(useCalibratedFiles) {
    //If we have CalibratedAnitaEvent then the constructor is just
    testEvent = new UsefulAnitaEvent(calEvent);
  }
  else {
    //If we have RawAnitaEvent then we have to specify the calibration option
    testEvent = new UsefulAnitaEvent(event,WaveCalType::kDefault,header);
  }
  int countH=0; int countV=0;
  int meanH=0,gotPhiH0=0;
  int meanV=0,gotPhiV0=0;
  for(int phi=0;phi<16;phi++) {
    std::cout << phi << "\t" << header->isInL3Pattern(phi,AnitaPol::kHorizontal) << "\t" << header->isInL3Pattern(phi,AnitaPol::kVertical) << "\n";
    if(header->isInL3Pattern(phi,AnitaPol::kHorizontal)) {
      countH++;
      meanH+=phi;
      if(phi==0 ) gotPhiH0=1;
      if(gotPhiH0 && phi>8) meanH-=16;
    }
      
     
    if(header->isInL3Pattern(phi,AnitaPol::kVertical)) {
      countV++;    
      meanV+=phi;
      if(phi==0 ) gotPhiV0=1;
      if(gotPhiV0 && phi>8) meanV-=16;
    }
  }
  int isHPolTest=0;
  int leftPhi=0;
  int centrePhi=0;
  int rightPhi=0;
  if(countH>countV) {
    std::cout << "H-Pol trigger\n";
    isHPolTest=1;
    meanH/=countH;
    if(meanH<0) meanH+=16;
    centrePhi=meanH;
    rightPhi=meanH+1;
    if(rightPhi>15) rightPhi-=16;
    leftPhi=meanH-1;
    if(leftPhi<0) leftPhi+=16;
    
    //    std::cout << meanH << "\t" << countH << "\n";
  }
  else if(countV>countH) {
    std::cout << "V-Pol trigger\n";    
    meanV/=countV;
    centrePhi=meanV;
    rightPhi=meanV+1;
    if(rightPhi>15) rightPhi-=16;
    leftPhi=meanV-1;
    if(leftPhi<0) leftPhi+=16;
    //    std::cout << meanV << "\t" << countV << "\n";
  }
  else {
    std::cout << "Confused Trigger\n";
    return ;
  }

  TGraph *grSeed[3][3];
  AnitaPol::AnitaPol_t testPol=AnitaPol::kVertical;
  if(isHPolTest) testPol=AnitaPol::kHorizontal;
  for(int ring=0;ring<3;ring++) {
    grSeed[ring][0]=testEvent->getGraph((AnitaRing::AnitaRing_t)ring,leftPhi,testPol);
    grSeed[ring][1]=testEvent->getGraph((AnitaRing::AnitaRing_t)ring,centrePhi,testPol);
    grSeed[ring][2]=testEvent->getGraph((AnitaRing::AnitaRing_t)ring,rightPhi,testPol);
  }




  
  Int_t numEnts=20000;
  Long64_t maxEntry=headTree->GetEntries();
  if(seedEntry+numEnts<maxEntry) maxEntry=seedEntry+numEnts;
  Long64_t startEntry=0;
  if(seedEntry-numEnts>0) startEntry=seedEntry-numEnts;

  Int_t starEvery=numEnts/100;
  if(starEvery==0) starEvery=1;  
  
  Long64_t countEvents=0;

  for(Long64_t entry=startEntry;entry<maxEntry;entry++) {
     if(entry%starEvery==0) std::cerr << "*";

     //Get header
     headTree->GetEntry(entry);
     eventNumber=header->eventNumber;
     {
       int countH=0; int countV=0;
       int meanH=0,gotPhiH0=0;
       int meanV=0,gotPhiV0=0;
       for(int phi=0;phi<16;phi++) {
	 if(header->isInL3Pattern(phi,AnitaPol::kHorizontal)) {
	   countH++;
	   meanH+=phi;
	   if(phi==0 ) gotPhiH0=1;
	   if(gotPhiH0 && phi>8) meanH-=16;
	 }
	 
	 
	 if(header->isInL3Pattern(phi,AnitaPol::kVertical)) {
	   countV++;    
	   meanV+=phi;
	   if(phi==0 ) gotPhiV0=1;
	   if(gotPhiV0 && phi>8) meanV-=16;
	 }
       }
       if(isHPolTest) {
	 if(countV>=countH) continue;
       }
       else {
	 if(countH>=countV) continue;
       }
       
       int leftPhi=0;
       int centrePhi=0;
       int rightPhi=0;
       if(isHPolTest) {
	 meanH/=countH;
	 if(meanH<0) meanH+=16;
	 centrePhi=meanH;
	 rightPhi=meanH+1;
	 if(rightPhi>15) rightPhi-=16;
	 leftPhi=meanH-1;
	 if(leftPhi<0) leftPhi+=16;
       }
       else {	 
	 meanV/=countV;
	 centrePhi=meanV;
	 rightPhi=meanV+1;
	 if(rightPhi>15) rightPhi-=16;
	 leftPhi=meanV-1;
	 if(leftPhi<0) leftPhi+=16;
       }
       
     }


     
     //Get event
     eventChain->GetEntry(entry);

     //Now we can make a UsefulAnitaEvent (or a UsefulAnitaEvent)
     UsefulAnitaEvent *realEvent=0;
     if(useCalibratedFiles) {
       //If we have CalibratedAnitaEvent then the constructor is just
       realEvent = new UsefulAnitaEvent(calEvent);
     }
     else {
       //If we have RawAnitaEvent then we have to specify the calibration option
       realEvent = new UsefulAnitaEvent(event,WaveCalType::kDefault,header);
     }


     TGraph *grTest[3][3];
     TGraph *grCor[3][3];
     AnitaPol::AnitaPol_t testPol=AnitaPol::kVertical;
     if(isHPolTest) testPol=AnitaPol::kHorizontal;
     for(int ring=0;ring<3;ring++) {
       grTest[ring][0]=realEvent->getGraph((AnitaRing::AnitaRing_t)ring,leftPhi,testPol);
       grTest[ring][1]=realEvent->getGraph((AnitaRing::AnitaRing_t)ring,centrePhi,testPol);
       grTest[ring][2]=realEvent->getGraph((AnitaRing::AnitaRing_t)ring,rightPhi,testPol);
     }
     


     sumPeakCorr=0;
     sumPeakCorrSq=0;

     for(int ring=0;ring<3;ring++) {
       for(int ant=0;ant<3;ant++) {
	 FFTtools::getWaveformSNR(grTest[ring][ant],peakToPeak[ring][ant],rms[ring][ant]);
	 grCor[ring][ant]=FFTtools::getCorrelationGraph(grSeed[ring][ant],grTest[ring][ant]);
	 Int_t peakBin=0;
	 FFTtools::getPeakSqVal(grCor[ring][ant],&peakBin);
	 peakCorr[ring][ant]=grCor[ring][ant]->GetY()[peakBin];
	 peakTime[ring][ant]=grCor[ring][ant]->GetX()[peakBin];
	 sumPeakCorr+=peakCorr[ring][ant];
	 sumPeakCorrSq+=peakCorr[ring][ant]*peakCorr[ring][ant];	
	 delete grTest[ring][ant];
	 delete grCor[ring][ant];
       }
     }
     for(int ant=0;ant<3;ant++) {
       peakToPeakBottomTop[ant]=peakToPeak[0][ant]/peakToPeak[2][ant];
     }
     
     if(realEvent) delete realEvent;
     corTree->Fill();
  }
  std::cerr << "\n";
  corTree->AutoSave();
  fpOut->Close();

 //  if(countEvents==0) {
//     std::cout << "No events written will tidy up after myself and delete "
// 	      << outName << "\n";
//     gSystem->Unlink(outName);
//   }
}
     
