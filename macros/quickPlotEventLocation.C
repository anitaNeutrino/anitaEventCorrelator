#include "AnitaConventions.h"
#include "UsefulAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "FFTtools.h"
#include <iostream>
#include <fstream>


const double TrueScaleLat=71;
const double CentralMeridian=0;
const double RadiusOfEarth=6378.1e3; //Metres
const double xOffest=375;
const double yOffset=312.5;
const double scale=271.5/2.19496e+06;
const double xSize=750;
const double ySize=625;


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


void quickPlotEventLocation(char *dirName, int run, int eventNum, double thetaWave, double phiWave) {
  RawAnitaHeader *header=0;
  Adu5Pat *pat=0;

  char headName[180];
  char gpsName[180];
  
  sprintf(headName,"%s/run%d/headFile%d.root",dirName,run,run);
  sprintf(gpsName,"%s/run%d/gpsFile%d.root",dirName,run,run);

  TFile *headFile = TFile::Open(headName);
  if(!headFile) {
    std::cerr << headName << "\n";
    return;
  }
  TFile *gpsFile = TFile::Open(gpsName);
  if(!gpsFile) {
    std::cerr << gpsName << "\n";
    return;
  }

  TTree *headTree = (TTree*) headFile->Get("headTree");
  if(!headTree) {
    std::cerr << "Couldn't open headTree\n";
    return;
  }
  headTree->SetBranchAddress("header",&header);


  TTree *adu5PatTree = (TTree*) gpsFile->Get("adu5PatTree");
  if(!adu5PatTree) {
    std::cerr << "Couldn't open adu5PatTree\n";
    return;
  }
  adu5PatTree->SetBranchAddress("pat",&pat);



  headTree->BuildIndex("eventNumber");
  Long64_t entry = headTree->GetEntryNumberWithIndex(eventNum);
  if(entry<0) {
    std::cerr << "Couldn't get event " << eventNum << "\n";
    return;
  }
  headTree->GetEntry(entry);

  adu5PatTree->BuildIndex("realTime");
  Long64_t adu5Entry = adu5PatTree->GetEntryNumberWithBestIndex(header->triggerTime);
  if(adu5Entry<0) {
    std::cerr << "Couldn't get GPS info\n";
    return;
  }
   
  adu5PatTree->GetEntry(adu5Entry);
  UsefulAdu5Pat usefulPat(pat);

  Double_t sourceLon,sourceLat;
  Double_t sourceX,sourceY;


  Double_t balloonLon,balloonLat;
  Double_t balloonX,balloonY;
  Double_t balloonAlt;

  std::cout << thetaWave << "\t" << phiWave << "\n";
  int solved=usefulPat.getSourceLonAndLatAltZero(phiWave,thetaWave,
						 sourceLon,sourceLat);
  if(!solved) {
    std::cerr << "Couldn't find location at altitude 0";
  }
  
  
  balloonLat=usefulPat.latitude;
  balloonLon=usefulPat.longitude;
  balloonAlt=usefulPat.altitude;
  getRelXYFromLatLong(balloonLat,balloonLon,balloonX,balloonY);

  getRelXYFromLatLong(sourceLat,sourceLon,sourceX,sourceY);
  std::cout << balloonLat << "\t" << balloonLon << "\t" << balloonAlt <<"\n";
  std::cout << sourceLat << "\t" << sourceLon << "\t" << 0 <<"\n";
  
   TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
   if(!canMap)
      canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
   canMap->Clear();
   canMap->SetLogz();
   canMap->SetTopMargin(0);
   canMap->SetBottomMargin(0);
   canMap->SetLeftMargin(0);
   canMap->SetRightMargin(0);
   TImage *img = TImage::Open("antarcticaIceMapBW.png");
   if (!img) {
      printf("Could not create an image... exit\n");
      return;
    }
   img->SetConstRatio(kFALSE);
   img->Draw("");
 
   TMarker *mark = new TMarker();//375./xSize,312.5./ySize,21);
   mark->SetMarkerColor(7);
   mark->SetMarkerStyle(23);
   mark->DrawMarker(balloonX,balloonY);
   mark->SetMarkerColor(6);
   mark->SetMarkerStyle(29);
   mark->DrawMarker(sourceX,sourceY);
   
 

}
