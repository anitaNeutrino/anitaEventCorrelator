#include "../../include/AnitaConventions.h"
#include "../../include/PrettyAnitaEvent.h"
#include "../../include/RawAnitaEvent.h"
#include "../../include/TimedAnitaHeader.h"
#include "../../include/PrettyAnitaHk.h"
#include "../../include/UsefulAdu5Pat.h"
#include "../../include/CorrelationSummary.h"
#include "../../include/AnitaGeomTool.h"
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



void correlationTreeLoop(int run);

  
void correlationTreeLoop(int run) {
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
   char eventName[FILENAME_MAX];
   char headerName[FILENAME_MAX];
   char hkName[FILENAME_MAX];
   char gpsName[FILENAME_MAX];
   char corrName[FILENAME_MAX];
   char outName[FILENAME_MAX];
   char baseDir[FILENAME_MAX];

   sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
 sprintf(eventName,"%s/run%d/eventFile%d.root",baseDir,run,run);
 sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
 sprintf(gpsName,"%s/run%d/gpsFile%d.root",baseDir,run,run);
 sprintf(corrName,"../../outFiles/corRun%d.root",run);
 sprintf(outName,"../../outFiles/deltaTFileSlowClock%d.root",run);

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
   
   TChain *eventChain = new TChain("eventTree");
   eventChain->Add(eventName);
   eventChain->SetBranchAddress("event",&event);
   
  //  TFile *fpHk = new TFile(hkName);
//    TTree *prettyHkTree = (TTree*) fpHk->Get("prettyHkTree");
//    prettyHkTree->SetBranchAddress("hk",&hk);
   
   TFile *fpHead = TFile::Open(headerName);
   TTree *headTree = (TTree*) fpHead->Get("headTree");
   headTree->SetBranchAddress("header",&header);
   
   TFile *fpGps = TFile::Open(gpsName);
   TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
   adu5PatTree->BuildIndex("realTime");
   adu5PatTree->SetBranchAddress("pat",&pat);

      TFile *fpCor = new TFile(corrName);
      TTree *corTree = (TTree*) fpCor->Get("corTree");
      corTree->SetBranchAddress("cor",&corSum);

   // TFile *fpGp = new TFile("/unix/anita/gpLogs/gp_log_all.root");
   //TTree *gpLogTree = (TTree*) fpGp->Get("gp_log");
   // Double_t gpTriggerTime;
   //Int_t antennaFlag, sgFreq;
   //gpLogTree->SetBranchAddress("trigger_time",&gpTriggerTime);
   //gpLogTree->SetBranchAddress("sg_freq",&sgFreq);
   //gpLogTree->SetBranchAddress("antenna",&antennaFlag);
   //gpLogTree->BuildIndex("trigger_time");

   Long64_t numEntries=headTree->GetEntries();
   int counter=0;

   TFile *fpOut = new TFile(outName,"RECREATE");
//    TH1F *histSimpleDtDiff  = new TH1F("histSimpleDtDiff","All DeltaT Diffs",2000,-10,10);
//    TH1F *histDtDiffSep[16][4][11];
//    char histName[180];
//    for(int ant=0;ant<16;ant++) {
//       for(int chip=0;chip<4;chip++) {
// 	 for(int cor=0;cor<11;cor++) {
// 	    sprintf(histName,"histDtDiffSep%d_%d_%d",ant,chip,cor);
// 	    histDtDiffSep[ant][chip][cor]= new TH1F(histName,histName,2000,-10,10);
// 	 }
//       }
//    }

   Long64_t entry=0;
   Int_t firstAnt,secondAnt,labChip,maxAnt;
   Double_t deltaT,deltaTExpected;
   Double_t phiWave, phiMaxAnt;
   Double_t corPeak, corRMS;
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
   deltaTTree->Branch("phiWave",&phiWave,"phiWave/D");
   deltaTTree->Branch("phiMaxAnt",&phiMaxAnt,"phiMaxAnt/D");
  

   Double_t thetaWave;

   for(entry=0;entry<numEntries;entry++) {
  
      //Stupidly most do this to be perfectly safe  
      headTree->GetEntry(entry);
     
      // if(header->triggerTimeNs*1e-9< 0.097 || header->triggerTimeNs*1e-9>0.1)
      if( (header->triggerTimeNs>0.3e6) || (header->triggerTimeNs<0.2e6) )  
      continue; 

      //Long64_t gpEntry=gpLogTree->GetEntryNumberWithBestIndex(header->triggerTime);
      //     cout << gpEntry << endl;
      //if(gpEntry>-1) 
      //	 gpLogTree->GetEntry(gpEntry);
      //else 
      //	 continue;
     
      //  if(antennaFlag!=1 || sgFreq!=-1) 
      //	 continue;

      //     cout << entry << "\t" << gpEntry << "\t" << header->triggerTime << "\t" << (UInt_t)gpTriggerTime << "\n";
     
      //  eventChain->GetEntry(entry);
      //prettyHkTree->GetEntry(entry);


     Long64_t bestEntry = adu5PatTree->GetEntryNumberWithBestIndex(header->triggerTime);
   if(bestEntry>-1) 
     adu5PatTree->GetEntry(bestEntry);
   else 
     continue;

   corTree->GetEntry(entry);

   // PrettyAnitaEvent realEvent(event,WaveCalType::kVTFullAGCrossCorClock,header);
      UsefulAdu5Pat usefulPat(pat);
      usefulPat.getThetaAndPhiWaveWillySeavey(thetaWave,phiWave);
      int ant=fGeomTool->getUpperAntNearestPhiWave(phiWave);

      //     cout << ant << "\t" << phiWave << "\t" << fGeomTool->getAntPhiPositionRelToAftFore(ant) << "\n";

      // labChip=realEvent.getLabChip(1);
      Double_t deltaTHere= 1. / (2.6*16.);
      //corSum = realEvent.getCorrelationSummary(ant,AnitaPol::kVertical,deltaTHere);

      //      usefulPat.getThetaAndPhiWaveWillySeavey(thetaWave,phiWave);
      //      int ant=fGeomTool->getUpperAntNearestPhiWave(phiWave);
      //      cout << ant << "\t" << corSum->firstAnt[1] << endl;
     
      for(int corInd=0;corInd<19;corInd++) {


	//replace taylor dome
	//	 deltaTExpected=usefulPat.getDeltaTWillySeavey(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	 deltaTExpected=usefulPat.getDeltaTTaylor(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);

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
	 deltaTTree->Fill();
	 //	cout << corSum->centreAntenna << "\t" << corSum->firstAnt[corInd] << "\t" << corSum->secondAnt[corInd] << "\t" << deltaTExpected << "\t" << corSum->maxCorTimes[corInd] << "\n";
      }


      //     cout << "Phi Compare:\t" << corSum->centreAntenna << "\t" << fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna)*TMath::RadToDeg() << "\t"
      //	  << usefulPat.getPhiWave()*TMath::RadToDeg() << endl;
      //     if(entry>10)
      //     cout << "Source:\t" << usefulPat.getSourceLongitude() << "\t" << usefulPat.getSourceLatitude() << "\t"
      //	  << usefulPat.getSourceAltitude() << "\n";
      //     cout << "Balloon:\t" << usefulPat.longitude << "\t" << usefulPat.latitude << "\t" << usefulPat.altitude << "\t" << usefulPat.heading << endl;
      counter++; 
      if(counter%100==0)
	 cerr << "*";
   }
//    histSimpleDtDiff->Write();
//    for(int ant=0;ant<16;ant++) {
//       for(int chip=0;chip<4;chip++) {
// 	 for(int cor=0;cor<11;cor++) {
// 	    histDtDiffSep[ant][chip][cor]->Write();
// 	 }
//       }
//    }
   deltaTTree->AutoSave();
   fpOut->Close();
   //   histSimpleDtDiff->Draw();
}

