#include "AnitaConventions.h"
#include "UsefulAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "CorrelationSummary.h"
#include "Adu5Pat.h"
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

const float TrueScaleLat=71;
const float CentralMeridian=0;
const float RadiusOfEarth=6378.1e3; //Metres
const float xOffest=375;
const float yOffset=312.5;
const float scale=271.5/2.19496e+06;
const float xSize=750;
const float ySize=625;

void getRelXYFromLatLong(double latitude, double longitude,
			 double &x, double &y)
{
    //Negative longitude is west
 //    //All latitudes assumed south
    double absLat=TMath::Abs(latitude);
    double r=RadiusOfEarth*TMath::Cos((90.-TrueScaleLat)*TMath::DegToRad())*TMath::Tan((90-absLat)*TMath::DegToRad());
    y=r*TMath::Cos(longitude*TMath::DegToRad());
    x=r*TMath::Sin(longitude*TMath::DegToRad());   

    y*=scale;
    y+=yOffset;
    y/=ySize;
    x*=scale;
    x+=xOffest;
    x/=xSize;
 
}

void loopOverCorTree(char *dirName, int run);

void loopOverCorTree()
{
   cout << "Usage: loopOverCorTree(run no.)\n";
   //  exampleLoopAllEvents(1030);
}
  

void loopOverCorTree(char *dirName, int run) {

   char corName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  sprintf(corName,"%s/corRun%d.root",dirName,run);
  sprintf(headerName,"/unix/anita1/newRootData/run%d/timedHeadFile%d.root",run,run);
  sprintf(gpsName,"/unix/anita1/newRootData/run%d/betterGpsFile%d.root",run,run);
  sprintf(outName,"sillyDST%d.root",run);

  TimedAnitaHeader *header =0;
  Adu5Pat *pat = 0;
  CorrelationSummary *corSum =0;
  

  TFile *fpHead = new TFile(headerName);
  TTree *headTree = (TTree*) fpHead->Get("headTree");
  headTree->SetBranchAddress("header",&header);

  TFile *fpPat = new TFile(gpsName);
  TTree *adu5PatTree = (TTree*) fpPat->Get("adu5PatTree");
  adu5PatTree->SetBranchAddress("pat",&pat);

  TFile *fpCor = new TFile(corName);
  TTree *corTree = (TTree*) fpCor->Get("corTree");
  corTree->SetBranchAddress("cor",&corSum);

  Double_t balloonLat, balloonLon, balloonAlt, heading;
  Double_t sourceLon,sourceLat,sourceX,sourceY,phiWave,thetaWave,phiWaveErr,thetaWaveErr;
  Double_t chiSq;
  Int_t eventNumber;
  UInt_t triggerTime,triggerTimeNs;

  TFile *fpOut = new TFile(outName,"RECREATE");
  TTree *sillyTree = new TTree("sillyTree","Tree of Reco Thingies");
  sillyTree->Branch("eventNumber",&eventNumber,"eventNumber/I");
  sillyTree->Branch("triggerTime",&triggerTime,"triggerTime/i");
  sillyTree->Branch("triggerTimeNs",&triggerTimeNs,"triggerTimeNs/i");
  sillyTree->Branch("balloonLat",&balloonLat,"balloonLat/D");
  sillyTree->Branch("balloonLon",&balloonLon,"balloonLon/D");
  sillyTree->Branch("balloonAlt",&balloonAlt,"balloonAlt/D");
  sillyTree->Branch("heading",&heading,"heading/D");
  sillyTree->Branch("phiWave",&phiWave,"phiWave/D");
  sillyTree->Branch("thetaWave",&thetaWave,"thetaWave/D");
  sillyTree->Branch("phiWaveErr",&phiWaveErr,"phiWaveErr/D");
  sillyTree->Branch("thetaWaveErr",&thetaWaveErr,"thetaWaveErr/D");
  sillyTree->Branch("sourceLat",&sourceLat,"sourceLat/D");
  sillyTree->Branch("sourceLon",&sourceLon,"sourceLon/D");
  sillyTree->Branch("sourceX",&sourceX,"sourceX/D");
  sillyTree->Branch("sourceY",&sourceY,"sourceY/D");
  sillyTree->Branch("chiSq",&chiSq,"chiSq/D");

  UInt_t numEntries=corTree->GetEntries();

  UInt_t count=0;
  UInt_t passCount=0;
  
  headTree->BuildIndex("eventNumber");
  
  cout << "There are " << numEntries << " in corTree of run " << run << endl;
  cout << "headTree->GetEntries()\t" << headTree->GetEntries() << endl;
  cout << "adu5PatTree->GetEntries()\t" << adu5PatTree->GetEntries() << endl;
  //  eventChain->GetEntry(0);
  for(UInt_t entry=0;entry<numEntries;entry++) {
    //Stupidly most do this to be perfectly safe  
     corTree->GetEntry(entry);

     Long64_t mainEntry=headTree->GetEntryNumberWithIndex(corSum->eventNumber);
     if(mainEntry<0) {
	std::cerr << "Couldn't get event : " << corSum->eventNumber << "\n";
	continue;
     }
     headTree->GetEntry(mainEntry);
     adu5PatTree->GetEntry(mainEntry);
     


     if(corSum->chiSq<5000  && corSum->thetaWave>=0) {
	UsefulAdu5Pat usefulPat(pat);
	int solved=usefulPat.getSourceLonAndLatAltZero(corSum->phiWave,corSum->thetaWave,
						       sourceLon,sourceLat);
	eventNumber=corSum->eventNumber;
	triggerTime=header->triggerTime;
	triggerTimeNs=header->triggerTimeNs;
	balloonLat=usefulPat.latitude;
	balloonLon=usefulPat.longitude;
	balloonAlt=usefulPat.altitude;
	heading=usefulPat.heading;
	phiWave=corSum->phiWave;
	phiWaveErr=corSum->phiWaveErr;
	thetaWave=corSum->thetaWave;
	thetaWaveErr=corSum->thetaWaveErr;
	getRelXYFromLatLong(sourceLat,sourceLon,sourceX,sourceY);
	chiSq=corSum->chiSq;
	if(solved) {
	   sillyTree->Fill();
	   passCount++;
	}



	//	cout << usefulPat.longitude << "\t" << usefulPat.latitude << "\t" << usefulPat.heading << "\t";
	//	cout << corSum->phiWave*TMath::RadToDeg() << "\t" << corSum->thetaWave*TMath::RadToDeg() << "\t";
	//	cout << sourceLon << "\t" << sourceLat << "\n";

     }

    count++;
  }
  sillyTree->AutoSave();
  fpOut->Close();
  cerr << endl;
  cout << "Processed " << passCount << " of " << count << " events.\n";
}

