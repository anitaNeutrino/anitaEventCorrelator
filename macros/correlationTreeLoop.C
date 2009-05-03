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
#include "TSystem.h"
#include <iostream>
#include <fstream>



void correlationTreeLoop(int run, char *baseDir, char *corTreeDir, char *outputDir);

  
void correlationTreeLoop(int run,char *baseDir, char *corTreeDir, char *outputDir) {
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
   //   fGeomTool->useKurtAnitaIINumbers(1);
   char eventName[FILENAME_MAX];
   char headerName[FILENAME_MAX];
   char gpsName[FILENAME_MAX];
   char corrName[FILENAME_MAX];
   char outName[FILENAME_MAX];

   //   sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
   sprintf(eventName,"%s/run%d/eventFile%d.root",baseDir,run,run);
   sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
   sprintf(gpsName,"%s/run%d/gpsEvent%d.root",baseDir,run,run);
   sprintf(corrName,"%s/corRun%d.root",corTreeDir,run);
   sprintf(outName,"%s/deltaTFile%d.root",outputDir,run);

   
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
   UInt_t expTaylorTime;
   TFile *fpCor = new TFile(corrName);
   TTree *corTree = (TTree*) fpCor->Get("corTree");
   corTree->SetBranchAddress("cor",&corSum);
   corTree->SetBranchAddress("labChip",&labChip);
   corTree->SetBranchAddress("expTaylorTime",&expTaylorTime);

   Long64_t numEntries=corTree->GetEntries();
   int counter=0;

   TFile *fpOut = new TFile(outName,"RECREATE");

   Long64_t entry=0;
   UInt_t eventNumber, triggerTime, triggerTimeNs;
   Int_t firstAnt,secondAnt,maxAnt,corInd;
   Double_t deltaT,deltaTExpected;
   Double_t phiWave, phiMaxAnt;
   Double_t thetaWave;
   Double_t corPeak, corRMS;
   Double_t balloonLat, balloonLon, balloonAlt;
   Double_t heading,pitch,roll;
   Double_t deltaZ;
   Double_t deltaR;
   Double_t meanPhiAntPair;
   Double_t deltaPhiAntPair;
   


   TTree *deltaTTree = new TTree("deltaTTree","Tree of Delta T's");
   deltaTTree->Branch("entry",&entry,"entry/L");
   deltaTTree->Branch("firstAnt",&firstAnt,"firstAnt/I");
   deltaTTree->Branch("secondAnt",&secondAnt,"secondAnt/I");
   deltaTTree->Branch("maxAnt",&maxAnt,"maxAnt/I");
   deltaTTree->Branch("labChip",&labChip,"labChip/I");
   deltaTTree->Branch("deltaT",&deltaT,"deltaT/D");
   deltaTTree->Branch("deltaTExpected",&deltaTExpected,"deltaTExpected/D");
   deltaTTree->Branch("corPeak",&corPeak,"corPeak/D");
   deltaTTree->Branch("corRMS",&corRMS,"corRMS/D");
   deltaTTree->Branch("phiMaxAnt",&phiMaxAnt,"phiMaxAnt/D");
   deltaTTree->Branch("phiWave",&phiWave,"phiWave/D");
   deltaTTree->Branch("thetaWave",&thetaWave,"thetaWave/D");
   deltaTTree->Branch("eventNumber",&eventNumber,"eventNumber/i");
   deltaTTree->Branch("triggerTime",&triggerTime,"triggerTime/i");
   deltaTTree->Branch("triggerTimeNs",&triggerTimeNs,"triggerTimeNs/i");
   deltaTTree->Branch("corInd",&corInd,"corInd/I");
   deltaTTree->Branch("balloonLat",&balloonLat,"balloonLat/D");
   deltaTTree->Branch("balloonLon",&balloonLon,"balloonLon/D");
   deltaTTree->Branch("balloonAlt",&balloonAlt,"balloonAlt/D");
   deltaTTree->Branch("heading",&heading,"heading/D");
   deltaTTree->Branch("pitch",&pitch,"pitch/D");
   deltaTTree->Branch("roll",&roll,"roll/D");
   deltaTTree->Branch("deltaZ",&deltaZ,"deltaZ/D");
   deltaTTree->Branch("deltaR",&deltaR,"deltaR/D");
   deltaTTree->Branch("meanPhiAntPair",&meanPhiAntPair,"meanPhiAntPair/D");
   deltaTTree->Branch("deltaPhiAntPair",&deltaPhiAntPair,"deltaPhiAntPair/D");
   deltaTTree->Branch("expTaylorTime",&expTaylorTime,"expTaylorTime/i");

  // Double_t thetaWave;

   for(entry=0;entry<numEntries;entry++) { 
  

      corTree->GetEntry(entry);
      Long64_t headEntry=headTree->GetEntryNumberWithIndex(corSum->eventNumber);
      if(headEntry<0) 
	continue;
      headTree->GetEntry(headEntry);
     
      // if(header->triggerTimeNs*1e-9< 0.097 || header->triggerTimeNs*1e-9>0.1)

      triggerTimeNs=header->triggerTimeNs;
      triggerTime=header->triggerTime;
      eventNumber=header->eventNumber;

      adu5PatTree->GetEntry(headEntry);
      
      
   // PrettyAnitaEvent realEvent(event,WaveCalType::kVTFullAGCrossCorClock,header);
      balloonLat=pat->latitude;
      balloonLon=pat->longitude;
      balloonAlt=pat->altitude;
      heading=pat->heading;
      pitch=pat->pitch;
      roll=pat->roll;
	//Simon numbers     
      pat->pitch=0.64;
      pat->roll=0.14;

      //Test heading offset
//       pat->heading+=0.24;
//       if(pat->heading>=360) pat->heading-=360;
//       if(pat->heading<0) pat->heading+=360;

     // 

      //      pat->pitch=0.5;
      //      pat->roll=-0.1;


      //Kurt numbers 
      //      pat->pitch=-0.29;
      //      pat->roll=0.89;
//       //Test heading offset
//       pat->heading-=0.16;
//       if(pat->heading>=360) pat->heading-=360;
//       if(pat->heading<0) pat->heading+=360;

      UsefulAdu5Pat usefulPat(pat);
     
      for(corInd=0;corInd<35;corInd++) {

	 firstAnt=corSum->firstAnt[corInd];
	 secondAnt=corSum->secondAnt[corInd];
	 deltaT=corSum->maxCorTimes[corInd];
	 maxAnt=corSum->centreAntenna;
	 phiMaxAnt=fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna)*TMath::RadToDeg();
	 
	 //Default values
	 deltaTExpected=usefulPat.getDeltaTWillySeavey(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	 deltaTExpected=usefulPat.getDeltaTTaylor(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	 

	 phiWave=usefulPat.getPhiWave();
	 thetaWave=usefulPat.getThetaWave();
	 corPeak=corSum->maxCorVals[corInd];
	 corRMS=corSum->rmsCorVals[corInd];
	 
	 meanPhiAntPair=fGeomTool->getMeanAntPairPhiRelToAftFore(firstAnt,secondAnt);
	 deltaPhiAntPair=fGeomTool->getPhiDiff(meanPhiAntPair,phiWave);
	 
	 //	 std::cout << phiWave << "\t" << meanPhiAntPair << "\t" << deltaPhiAntPair << "\n";

	 //Convert to degrees
	 deltaPhiAntPair*=TMath::RadToDeg();
	 phiWave*=TMath::RadToDeg();
	 thetaWave*=TMath::RadToDeg();	       

	 //Actually fill the tree
	 deltaTTree->Fill();
      }


      counter++; 
      if(counter%100==0)
	 cerr << "*";
   }
   deltaTTree->AutoSave();
   fpOut->Close();
   //   histSimpleDtDiff->Draw();
}

