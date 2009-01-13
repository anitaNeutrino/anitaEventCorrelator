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
#include "TAxis.h"

#define NUM_BINS_THETA 180
#define THETA_MAX 0
#define NUM_BINS_PHI 360
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

Int_t saturatedChannel[90];

void startCorrelation(int run,int eventNumber){

  //int entry=65000;

  char headName[FILENAME_MAX];
  char eventName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  RawAnitaEvent *evPtr = 0;
  RawAnitaHeader *hdPtr = 0;
  Adu5Pat *patPtr =0;

  sprintf(headName,"/TBdata/anita/antarctica08/webPlotter/events/root/run%d/headFile%d.root",run,run);
  sprintf(eventName,"/TBdata/anita/antarctica08/webPlotter/events/root/run%d/eventFile%d.root",run,run);
  sprintf(gpsName,"/TBdata/anita/antarctica08/webPlotter/events/root/run%d/gpsFile%d.root",run,run);

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
  headTree->BuildIndex("eventNumber");
  adu5PatTree->SetBranchAddress("pat",&patPtr);
  adu5PatTree->BuildIndex("realTime");

  int entry=headTree->GetEntryNumberWithBestIndex(eventNumber);

  eventTree->GetEntry(entry);
  headTree->GetEntry(entry);



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
  Int_t patEntry;
  patEntry = adu5PatTree->GetEntryNumberWithBestIndex(hdPtr->realTime);
  adu5PatTree->GetEntry(patEntry);

  std::cout << patPtr->realTime << " " << hdPtr->realTime << std::endl;

  if(patPtr->realTime != hdPtr->realTime){
    std:: cout << "pat time doesn't match head time, pat realTime: " << patPtr->realTime << " head realTime: " << hdPtr->realTime << std::endl;
    //return;
  }


  for(int chan=0;chan<90;chan++){
    if(evPtr->xMax[chan]>130 || evPtr->xMin[chan]<-130){
      saturatedChannel[chan]=1;
      std::cout << "saturated " << chan << std::endl;
    }
    else saturatedChannel[chan]=0;
  }


  //entry++;
  TH2D *crossCorrelation = crossCorrelate(evPtr,hdPtr,patPtr);
 
  sprintf(headName,"crossCorrCan");
  TCanvas *crossCorrCan = (TCanvas*)gROOT->FindObject(headName);
  if(!crossCorrCan)
    crossCorrCan = new TCanvas(headName,headName,800,400);
  crossCorrCan->Clear();
  //sumCrossCorrs->Draw("aitoff");
  crossCorrelation->Draw("colz");
  TAxis *phiSectors = new TAxis(16,1,16);
  phiSectors->Draw();

  //std::cout << "event " << hdPtr->eventNumber << " time " << 

  double theta,phi;

  getSignalDirection(crossCorrelation,phi,theta);

  std::cout << "phi " << phi << " theta " << theta << std::endl;

  plotAnitaEventMap(patPtr,phi,theta);

}



TGraph *getCorrelation(TGraph *gr1,TGraph *gr2){
  return FFTtools::getCorrelationGraph(gr1,gr2);
}


void getTriggeredPhi(RawAnitaHeader *hdPtr,int triggeredPhi[16]){
  int numPhiTrigs=0; 
  int extraPhi[16]={0}; 
  for(int phi=0;phi<16;phi++){
    if(hdPtr->l3TrigPattern & (1 << phi)){
      triggeredPhi[phi]=1;
      numPhiTrigs++;
    }
    else triggeredPhi[phi]=0;
  }
  if(numPhiTrigs==1){
    for(int phi=0;phi<16;phi++){
      if(triggeredPhi[phi]){
	if(phi+1<16)
	  triggeredPhi[phi+1]=1;
	else
	  triggeredPhi[0]=1;
	if(phi-1>0)
	  triggeredPhi[phi-1]=1;
	else
	  triggeredPhi[15]=1;
	break;
      }
    }
  }
  /*
  if(numPhiTrigs<3){
    for(int phi=0;phi<16;phi++){
      if(triggeredPhi[phi]){
	if(phi+1>15) extraPhi[0]=1;
	else extraPhi[phi+1]=1;
	if(phi-1<0) extraPhi[15]=1;
	else extraPhi[phi-1]=1;
      }
    }
    for(int phi=0;phi<16;phi++){
      if(extraPhi[phi]) triggeredPhi[phi]=1;
    }
  }
  */

}


void getTriggeredAnt(int triggeredPhi[16],int triggeredAnt[32]){
  int phi;
  int chanIndex;
  for(int ant=0;ant<32;ant++){
    phi=AnitaGeomTool::getPhiFromAnt(ant);
    if(triggeredPhi[phi]) triggeredAnt[ant]=1;
    else triggeredAnt[ant]=0;
    chanIndex=AnitaGeomTool::getChanIndexFromAntPol(ant,AnitaPol::kVertical);
    if(saturatedChannel[chanIndex]==1){
      triggeredAnt[ant]=0;
    }
    std::cout << "ant " << ant+1 << " phi " << phi+1 << " triggered? " << triggeredAnt[ant] << " " << triggeredPhi[phi] << std::endl;
  }
}


void setupCosSinArray(double thetaArray[NUM_BINS_THETA],double phiArray[NUM_BINS_PHI],double cosThetaArray[NUM_BINS_THETA],double sinThetaArray[NUM_BINS_THETA],double cosPhiArray[NUM_BINS_PHI],double sinPhiArray[NUM_BINS_PHI]){
  for(int i=0;i<NUM_BINS_PHI;i++){
    phiArray[i] = (i+0.5) * 2*PI/NUM_BINS_PHI;
    cosPhiArray[i] = cos(phiArray[i]);
    sinPhiArray[i] = sin(phiArray[i]);
  }
  for(int i=0;i<NUM_BINS_THETA;i++){
    thetaArray[i] = (i+0.5) * PI/NUM_BINS_THETA - PI/2;
    //thetaArray[i] = (i+0.5) * PI/(2*NUM_BINS_THETA) - PI/2;
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

  TGraph *grTemp[32]={0};
  TGraph *grTemp2[32]={0};
  TGraph *grTemp3[32]={0};
  TGraph *grInt[32]={0};
  TGraph *grCorr[528]={0};
  double firstCorrTime[528];
  double dummyY;
  Double_t deltaT=1/(2.6*8.);

  for(int ant=0;ant<32;ant++){
    if(triggeredAnt[ant]){
      //grInt[ant]=realEvent.getGraph(AnitaGeomTool::getChanIndexFromAntPol(ant,AnitaPol::kVertical));
      grTemp[ant]=realEvent.getGraph(AnitaGeomTool::getChanIndexFromAntPol(ant,AnitaPol::kVertical));
      grTemp2[ant]=FFTtools::simpleNotchFilter(grTemp[ant],370.,520.);
      grTemp3[ant]=FFTtools::simpleNotchFilter(grTemp2[ant],920.,980.);
      grInt[ant]=FFTtools::simpleNotchFilter(grTemp3[ant],1120,1220);
    }
  }
  int arrayRef=0;

  char canName[FILENAME_MAX];
  TCanvas *canny[528];

  for(int ant1=0;ant1<32;ant1++){
    for(int ant2=ant1;ant2<32;ant2++){

      if(triggeredAnt[ant1] && triggeredAnt[ant2] && ant1!=ant2){
	grCorr[arrayRef]=FFTtools::getInterpolatedCorrelationGraph(grInt[ant1],grInt[ant2],deltaT);
	grCorr[arrayRef]->GetPoint(0,firstCorrTime[arrayRef],dummyY);
	/*
	sprintf(canName,"can%d",arrayRef);
	canny[arrayRef] = new TCanvas(canName,canName,800,600);
	canny[arrayRef]->Divide(1,3);
	canny[arrayRef]->cd(1);
	sprintf(canName,"ant %d (phi %d)",ant1+1,AnitaGeomTool::getPhiFromAnt(ant1)+1);
	grInt[ant1]->SetTitle(canName);
	grInt[ant1]->Draw("al");
	canny[arrayRef]->cd(2);
	sprintf(canName,"ant %d (phi %d)",ant2+1,AnitaGeomTool::getPhiFromAnt(ant2)+1);
	grInt[ant2]->SetTitle(canName);
	grInt[ant2]->Draw("al");
	canny[arrayRef]->cd(3);
	sprintf(canName,"corr %d %d",ant1+1,ant2+1);
	grCorr[arrayRef]->SetTitle(canName);
	grCorr[arrayRef]->Draw("al");
	*/
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

	    //deltaTarray[phi][theta] = usefulPat.getDeltaTExpected(ant1,ant2,cosPhiArray[phi],sinPhiArray[phi],cosThetaArray[theta],sinThetaArray[theta]);
	    deltaTarray[phi][theta] = usefulPat.getDeltaTExpected(ant1,ant2,phiArray[phi],thetaArray[theta]);

	    getPoint=static_cast<int>((deltaTarray[phi][theta]-firstCorrTime[arrayRef])*1/deltaT);

	    grCorr[arrayRef]->GetPoint(getPoint,xVal1,pointVal1);
	    grCorr[arrayRef]->GetPoint(getPoint+1,xVal2,pointVal2);

	    //std::cout << "delta T" << deltaTarray[phi][theta] << " xVal " << xVal1 << " xVal2 " << xVal2 << std::endl;

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
  //H2D *sumCrossCorrs = new TH2D(histName,histName,NUM_BINS_PHI,2,17,NUM_BINS_THETA,-90,90);
  TH2D *sumCrossCorrs = new TH2D(histName,histName,NUM_BINS_PHI,0,360,NUM_BINS_THETA,-90,90);
  double phiscale = (15./360.);
  for(int phi=0;phi<NUM_BINS_PHI;phi++){
    for(int theta=0;theta<NUM_BINS_THETA;theta++){
      //sumCrossCorrs->Fill(phiArray[phi]*180/PI*phiscale+2.,thetaArray[theta]*180./PI,correlationArray[phi][theta]);
      if(thetaArray[theta]*180./PI < 50. && thetaArray[theta]*180./PI > -50.)// && phiArray[phi]*180./PI >270)
        sumCrossCorrs->Fill(phiArray[phi]*180./PI,thetaArray[theta]*180./PI,correlationArray[phi][theta]);

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

  cout << "max bin " << x << " " << y << " " << z << endl;

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

  cout << "sig phi " << phi << " theta " << theta << endl; 

}


void plotAnitaEventMap(Adu5Pat *patPtr,double phi,double theta){

  double sourceLon,sourceLat,headLon,headLat,phi10Lon,phi10Lat,phi6Lon,phi6Lat,phi14Lon,phi14Lat,actualLat,actualLon,actual2Lat,actual2Lon;
  float xEvent,yEvent,xAnita,yAnita,anitaLat,anitaLon,anitaAlt,xHead,yHead,x10,y10,x6,y6,x14,y14,yActual,xActual,yActual2,xActual2;

  anitaLat = patPtr->latitude;
  anitaLon = patPtr->longitude;
  anitaAlt = patPtr->altitude;

  UsefulAdu5Pat usefulPat(patPtr);
  std::cout << "source " << std::endl;
  //int sourceLoc = usefulPat.getSourceLonAndLatAltZero((180-phi)/180.*PI,(theta)/180.*PI,sourceLon,sourceLat);
  int sourceLoc = usefulPat.getSourceLonAndLatAltZero((phi)/180.*PI,(theta)/180.*PI,sourceLon,sourceLat);
  std::cout << std::endl << "heading " << std::endl;
  int headLoc = usefulPat.getSourceLonAndLatAltZero(0./180.*PI,10./180.*PI,headLon,headLat);
  std::cout << std::endl << "phi 10 " << std::endl;
  int headLoc10 = usefulPat.getSourceLonAndLatAltZero(180./180.*PI,10./180.*PI,phi10Lon,phi10Lat);
  std::cout << std::endl << "phi 14 " << std::endl;
  int headLoc14 = usefulPat.getSourceLonAndLatAltZero(270./180.*PI,10./180.*PI,phi14Lon,phi14Lat);
  std::cout << std::endl << "phi 6 " << std::endl;
  int headLoc6 = usefulPat.getSourceLonAndLatAltZero(90./180.*PI,10./180.*PI,phi6Lon,phi6Lat);
  std::cout << std::endl << "actual 14.5 " << std::endl;
  int actualLoc = usefulPat.getSourceLonAndLatAltZero(231.5/180.*PI,14.5/180.*PI,actualLon,actualLat);
  std::cout << std::endl << "actual 4.5 " << std::endl;
  int actualLoc2 = usefulPat.getSourceLonAndLatAltZero(231.5/180.*PI,7.5/180.*PI,actual2Lon,actual2Lat);
  //int sourceLoc = usefulPat.getSourceLonAndLatAltZero((phi)/180.*PI,(TMath::PiOver2()-theta)/180.*PI,sourceLon,sourceLat);
  TImage *map = TImage::Open("/home/anita/eventCorrelator/macros/antarcticaIceMap.png");

  std::cout << "sourceLoc " << sourceLoc << " phi " << phi << " theta " << theta << " lon " << sourceLon << " lat " << sourceLat << std::endl;
  gStyle->SetMarkerColor(kBlack);
  //gStyle->SetMarkerSize(2);
  gStyle->SetTextSize(0.02);
  TMarker *anitaPos = new TMarker(xAnita,yAnita,23);

  getRelXYFromLatLong(anitaLat,anitaLon,xAnita,yAnita);
  getRelXYFromLatLong(static_cast<float>(sourceLat),static_cast<float>(sourceLon),xEvent,yEvent);
  getRelXYFromLatLong(static_cast<float>(headLat),static_cast<float>(headLon),xHead,yHead);
  getRelXYFromLatLong(static_cast<float>(phi10Lat),static_cast<float>(phi10Lon),x10,y10);
  getRelXYFromLatLong(static_cast<float>(phi14Lat),static_cast<float>(phi14Lon),x14,y14);
  getRelXYFromLatLong(static_cast<float>(phi6Lat),static_cast<float>(phi6Lon),x6,y6);
  getRelXYFromLatLong(static_cast<float>(actualLat),static_cast<float>(actualLon),xActual,yActual);
  getRelXYFromLatLong(static_cast<float>(actual2Lat),static_cast<float>(actual2Lon),xActual2,yActual2);

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


  TMarker *headingPos = new TMarker(xHead,yHead,29);
  TMarker *heading14Pos = new TMarker(x14,y14,29);
  TMarker *heading10Pos = new TMarker(x10,y10,29);
  TMarker *heading6Pos = new TMarker(x6,y6,29);
  TMarker *actualPos = new TMarker(xActual,yActual,29);
  TMarker *actual2Pos = new TMarker(xActual2,yActual2,29);
  headingPos->SetMarkerColor(kRed);
  heading14Pos->SetMarkerColor(kGray);
  heading10Pos->SetMarkerColor(kBlack);
  heading6Pos->SetMarkerColor(kViolet);
  actualPos->SetMarkerColor(kYellow+3);
  actual2Pos->SetMarkerColor(kYellow+2);

  headingPos->Draw("");
  heading14Pos->Draw("");
  heading10Pos->Draw("");
  heading6Pos->Draw("");
  actualPos->Draw("");
  actual2Pos->Draw("");
  anitaPos->DrawMarker(xAnita,yAnita);
  

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

  TMarker *eventPos = new TMarker(xEvent,yEvent,29);

  eventPos->Draw("");

  sprintf(label,"ANITA location: lat %f; long %f; alt %f, x %f, y %f",anitaLat,anitaLon,anitaAlt,xAnita,yAnita);
  positionLabel = new TLatex();
  positionLabel->DrawLatex(0.05,0.97,label);
  sprintf(label,"Event location: lat %f; long %f, x %f, y %f",sourceLat,sourceLon,xEvent,yEvent);
  positionLabel->DrawLatex(0.05,0.94,label);

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


/*
void distanceCalculator(Adu5Pat *patPtr,double &distance,double &timeOfFlight){

  //balloon location
  double lat1=patPtr->latitude*(3.1415/180.);//degrees to radians
  double lon1=patPtr->longitude*(3.1415/180.);
  double height1=patPtr->altitude;//meters

  //taylor dome location
  double lat2=-77.8803*(3.1415/180.);
  double lon2=158.45925*(3.1415/180.);
  double height2=2260.;//meters
  
  double x1,y1,z1,x2,y2,z2;
  
  double lat,lon,h;
  
  for (int i=0;i<2;i++){
    if (i==0){
      lat=lat1;
      lon=lon1;
      h=height1;
    }
    if (i==1){
      lat=lat2;
      lon=lon2;
      h=height2;
    }
    double r=6378137.0;//radius of earth in meters
    double f=1/298.257223563;//flattening factor
    double C=pow(cos(lat)*cos(lat)+(1-f)*(1-f)*sin(lat)*sin(lat),-0.5);
    double Q=(1-f)*(1-f)*C;
    double x=(r*C+h)*cos(lat)*cos(lon);
    double y=(r*C+h)*cos(lat)*sin(lon);
    double z=(r*Q+h)*sin(lat);
    if (i==0){
      x1=x; 
      y1=y; 
      z1=z;
    }
    if (i==1){
      x2=x; 
      y2=y; 
      z2=z;
    }
  }

  distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  cout<<"Distance From Location 1 to Location 2 (km): "<<distance/1000.<<endl;
  timeOfFlight=distance/(0.299792458);//in nanoseconds
  cout<<"Time of Flight (us): "<<timeOfFlight/1000.<<endl;

}
*/
