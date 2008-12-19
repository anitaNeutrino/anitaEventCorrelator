//////////////////////////////////
///// CrossCorrelation.C
///// A macro that does some cross correlations and
///// interferometric map calculations ... basics of this have been
///// stolen from Ryan
///// Author: Matthew Mottram
//////////////////////////////////

#include <iostream>
#include <fstream>

//ANITA Includes
#include "PrettyAnitaEvent.h"
#include "UsefulAdu5Pat.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "RawAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "FFTtools.h"
#include "FFTWComplex.h"



//ROOT Includes
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TImage.h"
#include "TMarker.h"
#include "TLatex.h"

#define NUM_BINS_THETA 90
#define NUM_BINS_PHI 180
#define PI 3.14159265

//runs from palestine to use: 3847 - 3853, 3793, 3797, 3798

TGraph *getCorrelation(TGraph *gr1,TGraph *gr2);
void getTriggeredPhi(RawAnitaHeader *hdPtr,int triggeredPhi[16]);
TH2D *crossCorrelate(RawAnitaEvent *evPtr,RawAnitaHeader *hdPtr,Adu5Pat *patPtr);
void startCorrelation(int run,int entry);
void setupCosSinArray(double thetaArray[NUM_BINS_THETA],double phiArray[NUM_BINS_PHI],double cosThetaArray[NUM_BINS_THETA],double sinThetaArray[NUM_BINS_THETA],double cosPhiArray[NUM_BINS_PHI],double sinPhiArray[NUM_BINS_PHI]);
void getSignalDirection(TH2D *crossCorrelation,double &phi,double &theta);
void plotAnitaEventMap(Adu5Pat *patPtr,double phi,double theta);
void getRelXYFromLatLong(float latitude,float longitude,float &x,float &y);

//for the map plot
const float TrueScaleLat=71;
const float CentralMeridian=0;
const float RadiusOfEarth=6378.1e3; //Metres
const float xOffest=375;
const float yOffset=312.5;
const float scale=271.5/2.19496e+06;
const float xSize=750;
const float ySize=625;


void startCorrelation(int run,int entry){

  char headName[FILENAME_MAX];
  char eventName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  RawAnitaEvent *evPtr = 0;
  RawAnitaHeader *hdPtr = 0;
  Adu5Pat *patPtr =0;

  sprintf(headName,"/Users/Matt/WORK/ANITA/rootFiles/run%d/headFile%d.root",run,run);
  sprintf(eventName,"/Users/Matt/WORK/ANITA/rootFiles/run%d/eventFile%d.root",run,run);
  sprintf(gpsName,"/Users/Matt/WORK/ANITA/rootFiles/run%d/gpsFile%d.root",run,run);

  TFile *eventFile = new TFile(eventName);
  TFile *headFile = new TFile(headName);
  TFile *fpGps = new TFile(gpsName);

  if(!fpGps){
    std::cout << "no GPS file\n";
    return;
  }

  TTree *adu5PatTree = (TTree*) fpGps->Get("adu5PatTree");
  TTree *eventTree = (TTree*)eventFile->Get("eventTree");
  TTree *headTree = (TTree*)headFile->Get("headTree");

  eventTree->SetBranchAddress("event",&evPtr);
  headTree->SetBranchAddress("header",&hdPtr);
  adu5PatTree->SetBranchAddress("pat",&patPtr);

  eventTree->GetEntry(entry);
  headTree->GetEntry(entry);

  UInt_t lastEvent;

  while(hdPtr->eventNumber!=10900006){
    if(lastEvent==hdPtr->eventNumber){
      std::cout << "no more entries in run" << std::endl;
      return;
    }
    lastEvent=hdPtr->eventNumber;
    entry++;
    eventTree->GetEntry(entry);
    headTree->GetEntry(entry);
  }

  /*
  //for(int phi=0;phi<16;phi++){
    while(hdPtr->triggerTimeNs<349.99e6 || hdPtr->triggerTimeNs>350.005e6 || hdPtr->l3TrigPattern==0){
      if(lastEvent==hdPtr->eventNumber){
	std::cout << "no more entries in run" << std::endl;
	return;
      }
      std::cout << "opened entry " << entry << " (event " << hdPtr->eventNumber << ") with triggerTimeNs " << hdPtr->triggerTimeNs << std::endl;
      lastEvent=hdPtr->eventNumber;
      entry++;
      eventTree->GetEntry(entry);
      headTree->GetEntry(entry);
      adu5PatTree->GetEntry(entry);
    }
    //}
    */

  std::cout << "opened entry " << entry << " (event " << hdPtr->eventNumber << ") with triggerTimeNs " << hdPtr->triggerTimeNs << std::endl;

  //get the pat ptr that corresponds to the timing of the event
  adu5PatTree->GetEntry(0);
  Int_t firstPatTime = patPtr->realTime;
  Int_t patEntry = hdPtr->realTime - firstPatTime;
  adu5PatTree->GetEntry(patEntry);

  std::cout << patPtr->realTime << " " << hdPtr->realTime << std::endl;

  if(patPtr->realTime != hdPtr->realTime){
    std:: cout << "pat time doesn't match head time, pat realTime: " << patPtr->realTime << " head realTime: " << hdPtr->realTime << std::endl;
    return;
  }

  entry++;
  TH2D *crossCorrelation = crossCorrelate(evPtr,hdPtr,patPtr);
 
  sprintf(headName,"crossCorrCan");
  TCanvas *crossCorrCan = (TCanvas*)gROOT->FindObject(headName);
  if(!crossCorrCan)
    crossCorrCan = new TCanvas(headName,headName,800,400);
  crossCorrCan->Clear();
  //sumCrossCorrs->Draw("aitoff");
  crossCorrelation->Draw("colz");

  double theta,phi;

  getSignalDirection(crossCorrelation,phi,theta);

  //std::cout << "phi " << phi << " theta " << theta << std::endl;

  plotAnitaEventMap(patPtr,phi,theta);

}


TGraph *getCorrelation(TGraph *gr1,TGraph *gr2){
  return FFTtools::getCorrelationGraph(gr1,gr2);
}


void getTriggeredPhi(RawAnitaHeader *hdPtr,int triggeredPhi[16]){
  for(int phi=0;phi<16;phi++){
    if(hdPtr->l3TrigPattern & (1 << phi)) triggeredPhi[phi]=1;
    else triggeredPhi[phi]=0;
  }
}


void getTriggeredAnt(int triggeredPhi[16],int triggeredAnt[32]){
  int phi;
  for(int ant=0;ant<32;ant++){
    phi=AnitaGeomTool::getPhiFromAnt(ant);
    if(triggeredPhi[phi]) triggeredAnt[ant]=1;
    else triggeredAnt[ant]=0;
    std::cout << "ant " << ant << " phi " << phi << " triggered? " << triggeredAnt[ant] << " " << triggeredPhi[phi] << std::endl;
  }
}


void setupCosSinArray(double thetaArray[NUM_BINS_THETA],double phiArray[NUM_BINS_PHI],double cosThetaArray[NUM_BINS_THETA],double sinThetaArray[NUM_BINS_THETA],double cosPhiArray[NUM_BINS_PHI],double sinPhiArray[NUM_BINS_PHI]){
  for(int i=0;i<NUM_BINS_PHI;i++){
    phiArray[i] = (i+0.5) * 2*PI/NUM_BINS_PHI - PI;
    cosPhiArray[i] = cos(phiArray[i]);
    sinPhiArray[i] = sin(phiArray[i]);
  }
  for(int i=0;i<NUM_BINS_THETA;i++){
    thetaArray[i] = (i+0.5) * PI/NUM_BINS_THETA - PI/2;
    cosThetaArray[i] = cos(thetaArray[i]);
    sinThetaArray[i] = sin(thetaArray[i]);
  }
}


TH2D *crossCorrelate(RawAnitaEvent *evPtr,RawAnitaHeader *hdPtr,Adu5Pat *patPtr){

  gStyle->SetPalette(1);

  int triggeredPhi[16];
  int triggeredAnt[32];

  double thetaArray[NUM_BINS_THETA];
  double phiArray[NUM_BINS_PHI];
  double cosThetaArray[NUM_BINS_THETA];
  double sinThetaArray[NUM_BINS_THETA];
  double cosPhiArray[NUM_BINS_PHI];
  double sinPhiArray[NUM_BINS_PHI];

  PrettyAnitaEvent realEvent(evPtr,WaveCalType::kVTFullAGCrossCorClock,hdPtr);
  UsefulAdu5Pat usefulPat(patPtr);

  getTriggeredPhi(hdPtr,triggeredPhi);
  getTriggeredAnt(triggeredPhi,triggeredAnt);
  setupCosSinArray(thetaArray,phiArray,cosThetaArray,sinThetaArray,cosPhiArray,sinPhiArray);

  TGraph *grInt[32]={0};
  TGraph *grCorr[528]={0};
  double firstCorrTime[528];
  double dummyY;
  Double_t deltaT=1/(2.6*8);

  for(int ant=0;ant<32;ant++){
    if(triggeredAnt[ant]){
      grInt[ant]=realEvent.getGraph(AnitaGeomTool::getChanIndexFromAntPol(ant,AnitaPol::kVertical));
    }
  }
  int arrayRef=0;

  for(int ant1=0;ant1<32;ant1++){
    for(int ant2=ant1;ant2<32;ant2++){

      if(triggeredAnt[ant1] && triggeredAnt[ant2] && ant1!=ant2){
	grCorr[arrayRef]=FFTtools::getInterpolatedCorrelationGraph(grInt[ant1],grInt[ant2],deltaT);
	grCorr[arrayRef]->GetPoint(0,firstCorrTime[arrayRef],dummyY);
	arrayRef++;
      }
      else{
	if(ant1!=ant2) arrayRef++;
      }

    }
  }

  double deltaTarray[NUM_BINS_PHI][NUM_BINS_THETA];
  double correlationArray[NUM_BINS_PHI][NUM_BINS_THETA];

  for(int phi=0;phi<NUM_BINS_PHI;phi++){
    for(int theta=0;theta<NUM_BINS_THETA;theta++){
      deltaTarray[phi][theta]=0.;
      correlationArray[phi][theta]=0.;
    }
  }

  int getPoint;
  double xVal1;
  double xVal2;
  double weight1;
  double weight2;
  double pointVal1;
  double pointVal2;

  arrayRef=0;
  for(int ant1=0;ant1<32;ant1++){
    for(int ant2=ant1;ant2<32;ant2++){

      if(triggeredAnt[ant1] && triggeredAnt[ant2] && ant1!=ant2){// && ant1==2 && ant2==21){

	for(int phi=0;phi<NUM_BINS_PHI;phi++){
	  for(int theta=0;theta<NUM_BINS_THETA;theta++){

	    deltaTarray[phi][theta] = usefulPat.getDeltaTExpected(ant1,ant2,cosPhiArray[phi],sinPhiArray[phi],cosThetaArray[theta],sinThetaArray[theta]);
	    //deltaTarray[phi][theta] = usefulPat.getDeltaTExpected(ant1,ant2,phiArray[phi],thetaArray[theta]);

	    getPoint=static_cast<int>((deltaTarray[phi][theta]-firstCorrTime[arrayRef])*1/deltaT);

	    grCorr[arrayRef]->GetPoint(getPoint,xVal1,pointVal1);
	    grCorr[arrayRef]->GetPoint(getPoint+1,xVal2,pointVal2);

	    weight1 = 1 - fabs(deltaTarray[phi][theta]-xVal1)/(deltaT);
	    weight2 = 1 - fabs(deltaTarray[phi][theta]-xVal2)/(deltaT);

	    correlationArray[phi][theta]+=(pointVal1*weight1+pointVal2*weight2);

	  }//theta
	}//phi

      }//if

      if(ant1!=ant2) arrayRef++;

    }//ant2
  }//ant1

  char histName[FILENAME_MAX];
  sprintf(histName,"sumCrossCorrs");
  TH2D *sumCrossCorrs = new TH2D(histName,histName,NUM_BINS_PHI,-PI*180/PI,PI*180/PI,NUM_BINS_THETA,-PI/2*180/PI,PI/2*180/PI);
  
  for(int phi=0;phi<NUM_BINS_PHI;phi++){
    for(int theta=0;theta<NUM_BINS_THETA;theta++){
      sumCrossCorrs->Fill(phiArray[phi]*180/PI,thetaArray[theta]*180/PI,correlationArray[phi][theta]);
    }
  }

  //sprintf(histName,"crossCorrCan");
  //TCanvas *crossCorrCan = new TCanvas(histName,histName,800,400);
  //sumCrossCorrs->Draw("aitoff");
  //sumCrossCorrs->Draw("colz");

  return sumCrossCorrs;

}



void getSignalDirection(TH2D *crossCorrelation,double &phi,double &theta){

  int x,y,z;
  double phiUpCont,phiDownCont,thetaUpCont,thetaDownCont,phiMidCont,thetaMidCont,phiTotCont,thetaTotCont;
  crossCorrelation->GetMaximumBin(x,y,z);

  phiUpCont=crossCorrelation->GetXaxis()->GetBinCenter(x+1)*crossCorrelation->GetBinContent(x+1,y);
  phiDownCont=crossCorrelation->GetXaxis()->GetBinCenter(x-1)*crossCorrelation->GetBinContent(x-1,y);
  phiMidCont=crossCorrelation->GetXaxis()->GetBinCenter(x)*crossCorrelation->GetBinContent(x,y);

  thetaUpCont=crossCorrelation->GetYaxis()->GetBinCenter(y+1)*crossCorrelation->GetBinContent(x,y+1);
  thetaDownCont=crossCorrelation->GetYaxis()->GetBinCenter(y-1)*crossCorrelation->GetBinContent(x,y-1);
  thetaMidCont=crossCorrelation->GetYaxis()->GetBinCenter(y)*crossCorrelation->GetBinContent(x,y);

  phiTotCont=crossCorrelation->GetBinContent(x-1,y)+crossCorrelation->GetBinContent(x,y)+crossCorrelation->GetBinContent(x+1,y);
  thetaTotCont=crossCorrelation->GetBinContent(x,y-1)+crossCorrelation->GetBinContent(x,y)+crossCorrelation->GetBinContent(x,y+1);

  phi=(phiUpCont+phiMidCont+phiDownCont)/phiTotCont;
  theta=(thetaUpCont+thetaMidCont+thetaDownCont)/thetaTotCont;

}


void plotAnitaEventMap(Adu5Pat *patPtr,double phi,double theta){

  double sourceLon,sourceLat;
  float xEvent,yEvent,xAnita,yAnita,anitaLat,anitaLon,anitaAlt;

  anitaLat = patPtr->latitude;
  anitaLon = patPtr->longitude;
  anitaAlt = patPtr->altitude;

  UsefulAdu5Pat usefulPat(patPtr);
  int sourceLoc = usefulPat.getSourceLonAndLatAltZero(phi,theta,sourceLon,sourceLat);
  TImage *map = TImage::Open("/home/anita/eventCorrelator/macros/antarcticaIceMap.png");

  gStyle->SetMarkerColor(kBlack);
  gStyle->SetTextSize(0.02);
  TMarker *anitaPos = new TMarker(xAnita,yAnita,22);

  getRelXYFromLatLong(anitaLat,anitaLon,xAnita,yAnita);
  getRelXYFromLatLong(static_cast<float>(sourceLat),static_cast<float>(sourceLon),xEvent,yEvent);

  TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
  if(!canMap)
     canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
  canMap->Clear();
  canMap->SetLogz();
  canMap->SetTopMargin(0);
  canMap->SetBottomMargin(0);
  canMap->SetLeftMargin(0);
  canMap->SetRightMargin(0);

  map->Draw("");
  anitaPos->Draw("");

  TLatex *positionLabel=0;
  char label[FILENAME_MAX];

  if(sourceLoc==0){
    if(anitaAlt<0){
      sprintf(label,"Could not get event position, ANITA below 0 altitude!");
    }
    else if(theta>0){
      sprintf(label,"Pointing upwards!  Cannot locate source at ground position");
    }
    else{
      sprintf(label,"Unkown error, cannot position source at 0 altitude");
    }
    positionLabel = new TLatex();
    positionLabel->DrawLatex(0.05,0.95,label);
    return;
  }

  TMarker *eventPos = new TMarker(xAnita,yAnita,29);

  eventPos->Draw("");

  sprintf(label,"ANITA location: lat %f; long %f; alt %f",anitaLat,anitaLon,anitaAlt);
  positionLabel = new TLatex();
  positionLabel->DrawLatex(0.05,0.97,label);
  sprintf(label,"Event location: lat %f; long %f",sourceLat,sourceLon);
  positionLabel->DrawLatex(0.05,0.94,label);

  //std::cout << "anita lat " << anitaLat << " lon " << anitaLon << " alt " << anitaAlt << " x " << xAnita << " y " << yAnita << std::endl;
  //std::cout << "source lat " << sourceLat << " lon " << sourceLon << " x " << xEvent << " y " << yEvent << std::endl;

}



void getRelXYFromLatLong(float latitude, float longitude,float &x, float &y)
{
    //Negative longitude is west
 //    //All latitudes assumed south
    float absLat=TMath::Abs(latitude);
    float r=RadiusOfEarth*TMath::Cos((90.-TrueScaleLat)*TMath::DegToRad())*TMath::Tan((90-absLat)*TMath::DegToRad());
    y=r*TMath::Cos(longitude*TMath::DegToRad());
    x=r*TMath::Sin(longitude*TMath::DegToRad());   

    y*=scale;
    y+=yOffset;
    y/=ySize;
    x*=scale;
    x+=xOffest;
    x/=xSize;
 
}
