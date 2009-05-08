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
#include "TMultiGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>


TCanvas *canSurf= new TCanvas("canSurf","canSurf");

double newOpt(double *par);
vector<double> leastSquares(double *xIn, double *yIn, int count);
  
double newOpt(double *par){

  int upperAntNums[NUM_PHI]={8,0,9,1,10,2,11,3,12,4,13,5,14,6,15,7};
  int lowerAntNums[NUM_PHI]={16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
  int lowerAntFromUpper[NUM_PHI]={17,19,21,23,25,27,29,31,16,18,20,22,24,26,28,30};
  int nadirAntNums[NUM_PHI]={32,-1,33,-1,34,-1,35,-1,36,-1,37,-1,38,-1,39,-1};
  int bottomFromNadir[8]={16,18,20,22,24,26,28,30};

  Double_t deltaR[40]={0};  
  Double_t deltaZ[40]={0};  
  Double_t deltaPhi[40]={0};  
  Double_t deltaHeading[1]={0};
  double  deltaTArrayMod[40]={0};
 
  int count3 = 0;
  int count17 = 0;

  double deltaTArray[60000] = {0};
  double phiAngleArray[60000]= {0};

  double deltaTArrayUpper[60000] = {0};
  double phiAngleArrayUpper[60000]= {0};
  double thetaAngleArrayUpper[60000]= {0};

  double deltaTArrayLoopUpper[60000] ={0};
  double phiAngleArrayLoopUpper[60000] = {0};
  double thetaAngleArrayLoopUpper[60000] = {0};

  double deltaTArrayLoop[60000] ={0};
  double phiAngleArrayLoop[60000] = {0};

  int middleAnt; 
  int leftAnt,rightAnt;

  int countArray = 0;
  int countArrayUpper = 0;


  for(int i = 0; i<32;i++){
    deltaR[i]=par[i];   
    deltaZ[i]=par[i+32];
    deltaPhi[i]=par[i+64];
    deltaTArrayMod[i]=par[i+96];
    cout << i << "  " <<deltaTArrayMod[i] << " " << deltaR[i] << " " << deltaPhi[i] << " " << deltaZ[i] << endl;

  }


  double theReturn = 0;
  double sumMean = 0;
  double sumMean2 = 0;
  double sumMean3 = 0;
  double sumMean4 = 0;
  double sumMean5 = 0;
  double sumMean6 = 0;
  int count8 = 0;
  int count18 = 0;
  double sumGrads = 0;

  TMultiGraph *myMG = new TMultiGraph;
  TMultiGraph *myMG3 = new TMultiGraph;
  TMultiGraph *myMG2 = new TMultiGraph;;

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corrName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  char baseDir[FILENAME_MAX];
  char *corTreeDir = "../../../Outfiles";
  double dummyArray[40] ={0}; 
  TGraph *tempAntGraph;
  TGraph *tempAntGraph2;
  TGraph *tempAntGraph3;

 
 vector<vector<double> > phiAngle;
  vector<vector<double> > deltaTVec;
  vector<vector<double> > thetaAngle;
  vector<vector<int> > firstAntVec;
  vector<vector<int> > secondAntVec;

  vector<double> phiAngleArray2;
  vector<double> deltaTArray2;

  vector<double> phiAngleArray2Upper;
  vector<double> thetaAngleArray2Upper;
  vector<double> deltaTArray2Upper;

  vector<double> temp;
  vector<int> temp2;
  temp.push_back(0);
  temp2.push_back(0);



  int leftOpt, rightOpt;


  double meanPhi[40] = {0}; 
  meanPhi[0] =22.5-12.5;
  meanPhi[1] =67.5-12.5;
  meanPhi[2] =112.5-12.5;
  meanPhi[3] =157.5-12.5;
  meanPhi[4] =202.5-12.5;
  meanPhi[5] =247.5-12.5;
  meanPhi[6] =292.5-12.5;
  meanPhi[7] =337.5-12.5;
  meanPhi[8] =45-12.5;
  meanPhi[9] =90-12.5;
  meanPhi[10] =135-12.5;
  meanPhi[11] =180-12.5;
  meanPhi[12] =225-12.5;
  meanPhi[13] =270-12.5;
  meanPhi[14] =315-12.5;
  meanPhi[15] =360-12.5;

  meanPhi[16] = 22.5-12.5;
  meanPhi[17] = 45-12.5;
  meanPhi[18] = 67.5-12.5;
  meanPhi[19] = 90-12.5;
  meanPhi[20] = 112.5-12.5;
  meanPhi[21] = 135-12.5;
  meanPhi[22] = 157.5-12.5;
  meanPhi[23] = 180-12.5;
  meanPhi[24] = 202.5-12.5;
  meanPhi[25] = 225-12.5;
  meanPhi[26] = 247.5-12.5;
  meanPhi[27] = 270-12.5;
  meanPhi[28] = 292.5-12.5;
  meanPhi[29] = 315-12.5;
  meanPhi[30] = 337.5-12.5;
  meanPhi[31] = 360-12.5;


  meanPhi[32] = 22.5-12.5;
  meanPhi[33] = 67.5-12.5;
  meanPhi[34] = 112.5-12.5;
  meanPhi[35] = 157.5-12.5;
  meanPhi[36] = 202.5-12.5;
  meanPhi[37] = 247.5-12.5;
  meanPhi[38] = 292.5-12.5;
  meanPhi[39] = 337.5-12.5;

//   for(int i =0; i < 32; i++){
//     phiAngleArray2.push_back(0);
//     deltaTArray2.push_back(0);
//     phiAngleArray2Upper.push_back(temp);
//     deltaTArray2Upper.push_back(temp);
//     thetaAngleArray2Upper.push_back(temp);

//   }

  for(int loop = 1; loop <4; loop++){
    int run = 16+loop;
    //int run = 18;

    canSurf->cd(loop+1);

    phiAngle.clear();
    thetaAngle.clear();
    deltaTVec.clear();
    firstAntVec.clear();
    secondAntVec.clear();

    for(int i = 0; i < 32; i++){
      phiAngle.push_back(temp);
      thetaAngle.push_back(temp);
      deltaTVec.push_back(temp);
      firstAntVec.push_back(temp2);
      secondAntVec.push_back(temp2);
    }

    //sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
    sprintf(baseDir,"/Users/simonbevan/Desktop/");
    sprintf(eventName,"%s/run%d/eventFile%d.root",baseDir,run,run);
    sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
    sprintf(gpsName,"%s/run%d/gpsFile%d.root",baseDir,run,run);
    sprintf(corrName,"%s/corRun%d.root",corTreeDir,run);

   
    RawAnitaEvent *event = 0;
    PrettyAnitaHk *hk = 0;
   
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
    TFile *fpCor = new TFile(corrName);
    TTree *corTree = (TTree*) fpCor->Get("corTree");
    corTree->SetBranchAddress("cor",&corSum);
    corTree->SetBranchAddress("labChip",&labChip);

    Long64_t numEntries=corTree->GetEntries();
    int counter=0;

    Long64_t entry=0;
    UInt_t eventNumber, triggerTime, triggerTimeNs;
    Int_t firstAnt,secondAnt,maxAnt,corInd;
    Double_t deltaT,deltaTExpected;
    Double_t phiWave, phiMaxAnt;
    Double_t corPeak, corRMS;
    Double_t balloonLat, balloonLon, balloonAlt;
    Double_t heading,pitch,roll;
  

    Double_t thetaWave;

    for(entry=0;entry<numEntries;entry++) {

      corTree->GetEntry(entry);
      Long64_t headEntry=headTree->GetEntryNumberWithIndex(corSum->eventNumber);
      if(headEntry<0) 
	continue;
      headTree->GetEntry(headEntry);
     
      if( (header->triggerTimeNs>0.5e6) || (header->triggerTimeNs<0.2e6) )  
	continue; 

      triggerTimeNs=header->triggerTimeNs;
      triggerTime=header->triggerTime;
      eventNumber=header->eventNumber;
  

      Long64_t bestEntry = adu5PatTree->GetEntryNumberWithBestIndex(header->triggerTime);
      if(bestEntry>-1) 
	adu5PatTree->GetEntry(bestEntry);
      else 
	continue;
      
      balloonLat=pat->latitude;
      balloonLon=pat->longitude;
      balloonAlt=pat->altitude;
      heading=pat->heading;
      pat->pitch=0.64;
      pat->roll=0.14;

      pitch=pat->pitch;
      roll=pat->roll;

      UsefulAdu5Pat usefulPat(pat);
      
      for(corInd=0;corInd<19;corInd++) {
	
	firstAnt=corSum->firstAnt[corInd];
	secondAnt=corSum->secondAnt[corInd];

	//replace taylor dome
	
	usefulPat.fSourceLongitude=0;
	//	deltaTExpected=usefulPat.getDeltaTTaylor(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	deltaTExpected=usefulPat.getDeltaTTaylorOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd],deltaR,deltaZ,deltaPhi);
	
	deltaT=corSum->maxCorTimes[corInd];
	maxAnt=corSum->centreAntenna;
	phiWave=usefulPat.getPhiWave()*TMath::RadToDeg();
	thetaWave=usefulPat.getThetaWave()*TMath::RadToDeg();
	phiMaxAnt=fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna)*TMath::RadToDeg();
	corPeak=corSum->maxCorVals[corInd];
	corRMS=corSum->rmsCorVals[corInd];

	Double_t meanPhiAntPair=fGeomTool->getMeanAntPairPhiRelToAftFore(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	Double_t deltaPhiAntPair=fGeomTool->getPhiDiff(meanPhiAntPair,usefulPat.getPhiWave());

	if((deltaT - deltaTExpected)*(deltaT - deltaTExpected) < 1 && (corPeak/corRMS)>8 && TMath::Abs(deltaPhiAntPair)*TMath::RadToDeg()<20){
	  phiAngle[0].push_back(phiWave);
	  thetaAngle[0].push_back(thetaWave);
	  deltaTVec[0].push_back(deltaT - deltaTExpected + deltaTArrayMod[firstAnt] - deltaTArrayMod[secondAnt]);
	  firstAntVec[0].push_back(firstAnt);
	  secondAntVec[0].push_back(secondAnt);  
	}
      }

      counter++; 
    
    }

   delete event;
       delete hk; 
       delete header;
       delete pat;
       delete corSum;

       delete fpHead;
       delete fpGps ;
       delete fpCor;
 
    //fill arrays
    //for(int ants = par[0]; ants < par[0]+1; ants++){

  sumMean = 0;
      sumMean2 = 0;
      sumMean3 = 0;
      sumMean4 = 0;
      sumMean5 = 0;
      sumMean6 = 0;
      sumGrads = 0;

    for(int ants = 0; ants < 32; ants++){

 
      int count = 0;
      int count2 = 0;
      int count16 = 0;
      count17 = 0;
      count3 = 0;
     
      if(ants <32){
	fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
      }else{
	leftAnt = ants;
	rightAnt = ants +1;
	if(ants == 39){
	  leftAnt = ants;
	  rightAnt = 32;
	}
      }    

      for(int events = 1; events < phiAngle[0].size(); events++){
	int firstAntTemp = (int)firstAntVec[0][events];
	int secondAntTemp = (int)secondAntVec[0][events];
	int rightTemp = int(rightAnt);

	int aboveTemp = 0;

	if(ants <16){
	  aboveTemp = lowerAntFromUpper[ants];
	}else{
	  aboveTemp = upperAntNums[ants-16];
	}
  
	

	if(firstAntTemp < 32){
      
 	  if( ((firstAntTemp == ants) &&  (secondAntTemp == rightTemp))){
	    //if((firstAntTemp == ants) &&  (secondAntTemp == rightTemp)){
	      deltaTArray[count] = deltaTVec[0][events];
	      phiAngleArray[count] = phiAngle[0][events];
	      count++;

	  } else if(((firstAntTemp == ants) &&  (secondAntTemp == aboveTemp))){
	    // 	    //if((firstAntTemp == ants) &&  (secondAntTemp == rightTemp)){

	  	      deltaTArrayUpper[count16] = deltaTVec[0][events];
	      phiAngleArrayUpper[count16] = phiAngle[0][events];
	      thetaAngleArrayUpper[count16] = thetaAngle[0][events];

	      count16++;   
	    
	  }
 
 	}
    
    
      }
      countArray = count;
      countArrayUpper = count16;
  
   
      count3 = 0;
      if(ants <32){
	fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
      }else{
	leftAnt = ants -1;
	rightAnt = ants +1;
	if(ants == 39){
	  leftAnt = ants - 1;
	  rightAnt = 32;
	}
      }    

    


      for(int events = 0; events < countArray; events++){
	  	    
	phiAngleArray2.push_back(phiAngleArray[events]);
 	deltaTArray2.push_back(deltaTArray[events]);
	count3++;

      }
	
      for(int events = 0; events < countArrayUpper; events++){
	  
	thetaAngleArray2Upper.push_back(thetaAngleArrayUpper[events]);
	phiAngleArray2Upper.push_back(phiAngleArrayUpper[events]);
       	deltaTArray2Upper.push_back(deltaTArrayUpper[events]);	
	count17++;
	
      }
	
      count8 = 0;
      count18 = 0;

      for(int events = 1; events < phiAngleArray2.size(); events++){

	if( deltaTArrayLoop[count8]<1){  
	  deltaTArrayLoop[count8] = deltaTArray2[events];
	  phiAngleArrayLoop[count8] = phiAngleArray2[events];
	  count8++;
      
	}
      
      }

      for(int events = 1; events < phiAngleArray2Upper.size(); events++){

	if( deltaTArrayLoopUpper[count18]<1){      
	  deltaTArrayLoopUpper[count18] = deltaTArray2Upper[events];
	  phiAngleArrayLoopUpper[count18] = phiAngleArray2Upper[events];
	  thetaAngleArrayLoopUpper[count18] = thetaAngleArray2Upper[events];
	  count18++;
      
	}
      
      }


      // cout << ants << count8 << "  " << count18 << endl;
	
      //   cout << count8 << endl;
       

      if(count8==0){

	tempAntGraph  = new TGraph(1, dummyArray, dummyArray);
	//tempAntGraph2  = new TGraph(1, dummyArray, dummyArray);
	//tempAntGraph3  = new TGraph(1, dummyArray, dummyArray);
      }else{
	tempAntGraph  = new TGraph(count8-1,  phiAngleArrayLoop, deltaTArrayLoop);
	//	tempAntGraph2  = new TGraph(count18-1,  phiAngleArrayLoopUpper, deltaTArrayLoopUpper);      
	//tempAntGraph3  = new TGraph(count18-1,  thetaAngleArrayLoopUpper, deltaTArrayLoopUpper);      

	if(ants < 32){  
	
	  //	canSurf->cd(1);
       

	  tempAntGraph->SetMinimum(-0.5);
	  tempAntGraph->SetMaximum(0.5);
	  tempAntGraph->Draw("ap");    

	  tempAntGraph->SetMarkerStyle(1);
	  tempAntGraph->GetXaxis()->SetLimits(0,360);
    
	

	  if(ants < 16){
	    sumMean = sumMean + tempAntGraph->GetMean(2);
	    sumMean5 = sumMean5 + tempAntGraph->GetMean(2)*tempAntGraph->GetMean(2);
	  }else{

	    sumMean4 = sumMean4 + tempAntGraph->GetMean(2);
	    sumMean6 = sumMean6 + tempAntGraph->GetMean(2)*tempAntGraph->GetMean(2);
	  }



	  if(ants == 8 || ants == 16 || ants == 12 || ants == 24){
	    tempAntGraph->SetMarkerColor(8);
	  }

	  if(ants == 3 || ants == 7 || ants == 23 || ants == 31){
	    tempAntGraph->SetMarkerColor(1);
	  }
    
	  if(ants == 9 || ants == 13 || ants == 18 || ants == 26){
	    tempAntGraph->SetMarkerColor(2);
	  }
    
    
	  if(ants == 10 || ants == 14 || ants == 20 || ants == 28){
	    tempAntGraph->SetMarkerColor(3);
	  }
    
	  if(ants == 11 || ants == 15 || ants == 22 || ants == 30){
	    tempAntGraph->SetMarkerColor(4);
	  }
    
	  if(ants == 0 || ants == 4 || ants == 17 || ants == 25){
	    tempAntGraph->SetMarkerColor(5);
	  }


	  if(ants == 2 || ants == 6 || ants == 21 || ants == 29){
	    tempAntGraph->SetMarkerColor(6);
	  }

	  if(ants == 1 || ants == 5 || ants == 19 || ants == 27){
	    tempAntGraph->SetMarkerColor(7);
	  }
    
	  tempAntGraph->GetXaxis()->SetTitle("phi (degrees)");
	  tempAntGraph->GetYaxis()->SetTitle("actual - expected time");
	  myMG2->Add(tempAntGraph);
	  myMG2->Draw("p");

	  vector<double> myFit = leastSquares(phiAngleArrayLoop, deltaTArrayLoop, count8-1);
	  double slope = myFit[0];
	  double intercept = myFit[1];

	  double tempX[2] = {slope*(meanPhi[ants]-20)+intercept,slope*(meanPhi[ants]+20)+intercept};
	  double tempY[2] = {(meanPhi[ants]-20),(meanPhi[ants]+20)};
	  sumGrads = sumGrads + myFit[0]*myFit[0]*10000;


	}

	if(count18==0){
	  
	  //tempAntGraph  = new TGraph(1, dummyArray, dummyArray);
	tempAntGraph2  = new TGraph(1, dummyArray, dummyArray);
	tempAntGraph3  = new TGraph(1, dummyArray, dummyArray);
      }else{
	  //tempAntGraph  = new TGraph(count8-1,  phiAngleArrayLoop, deltaTArrayLoop);
	tempAntGraph2  = new TGraph(count18-1,  phiAngleArrayLoopUpper, deltaTArrayLoopUpper);      
	tempAntGraph3  = new TGraph(count18-1,  thetaAngleArrayLoopUpper, deltaTArrayLoopUpper);      
	

	  tempAntGraph2->SetMinimum(-0.5);
	  tempAntGraph2->SetMaximum(0.5);
	  tempAntGraph2->Draw("ap");    

	  tempAntGraph2->SetMarkerStyle(1);
	  tempAntGraph2->GetXaxis()->SetLimits(0,360);
    
	  sumMean2 = sumMean2 + tempAntGraph2->GetMean(2)*tempAntGraph2->GetMean(2);

	  if(ants == 8 || ants == 16 || ants == 12 || ants == 24){
	    tempAntGraph2->SetMarkerColor(8);
	  }

	  if(ants == 3 || ants == 7 || ants == 23 || ants == 31){
	    tempAntGraph2->SetMarkerColor(1);
	  }
    
	  if(ants == 9 || ants == 13 || ants == 18 || ants == 26){
	    tempAntGraph2->SetMarkerColor(2);
	  }
    
    
	  if(ants == 10 || ants == 14 || ants == 20 || ants == 28){
	    tempAntGraph2->SetMarkerColor(3);
	  }
    
	  if(ants == 11 || ants == 15 || ants == 22 || ants == 30){
	    tempAntGraph2->SetMarkerColor(4);
	  }
    
	  if(ants == 0 || ants == 4 || ants == 17 || ants == 25){
	    tempAntGraph2->SetMarkerColor(5);
	  }


	  if(ants == 2 || ants == 6 || ants == 21 || ants == 29){
	    tempAntGraph2->SetMarkerColor(6);
	  }

	  if(ants == 1 || ants == 5 || ants == 19 || ants == 27){
	    tempAntGraph2->SetMarkerColor(7);
	  }
    
	  tempAntGraph2->GetXaxis()->SetTitle("phi (degrees)");
	  tempAntGraph2->GetYaxis()->SetTitle("actual - expected time");
	  myMG2->Add(tempAntGraph2);
	  myMG2->Draw("p");


	  tempAntGraph3->SetMinimum(-0.5);
	  tempAntGraph3->SetMaximum(0.5);
	  tempAntGraph3->Draw("ap");    

	  tempAntGraph3->SetMarkerStyle(1);
	  tempAntGraph3->GetXaxis()->SetLimits(0,360);
    
	  sumMean3 = sumMean3 + tempAntGraph3->GetMean(2)*tempAntGraph3->GetMean(2);

	  if(ants == 8 || ants == 16 || ants == 12 || ants == 24){
	    tempAntGraph3->SetMarkerColor(8);
	  }

	  if(ants == 3 || ants == 7 || ants == 23 || ants == 31){
	    tempAntGraph3->SetMarkerColor(1);
	  }
    
	  if(ants == 9 || ants == 13 || ants == 18 || ants == 26){
	    tempAntGraph3->SetMarkerColor(2);
	  }
    
    
	  if(ants == 10 || ants == 14 || ants == 20 || ants == 28){
	    tempAntGraph3->SetMarkerColor(3);
	  }
    
	  if(ants == 11 || ants == 15 || ants == 22 || ants == 30){
	    tempAntGraph3->SetMarkerColor(4);
	  }
    
	  if(ants == 0 || ants == 4 || ants == 17 || ants == 25){
	    tempAntGraph3->SetMarkerColor(5);
	  }


	  if(ants == 2 || ants == 6 || ants == 21 || ants == 29){
	    tempAntGraph3->SetMarkerColor(6);
	  }

	  if(ants == 1 || ants == 5 || ants == 19 || ants == 27){
	    tempAntGraph3->SetMarkerColor(7);
	  }
    
	  tempAntGraph3->GetXaxis()->SetTitle("phi (degrees)");
	  tempAntGraph3->GetYaxis()->SetTitle("actual - expected time");
	  myMG2->Add(tempAntGraph3);
	  myMG2->Draw("p");







	  vector<double> myFit = leastSquares(phiAngleArrayLoopUpper, deltaTArrayLoopUpper, count18-1);
	  double slope = myFit[0];
	  double intercept = myFit[1];	
       
	  	sumGrads = sumGrads + myFit[0]*myFit[0]*10000;


	  //cout << slope << endl;
      
	}
      } 
    

      phiAngleArray2.clear();
      deltaTArray2.clear();
      phiAngleArray2Upper.clear();
      thetaAngleArray2Upper.clear();
      deltaTArray2Upper.clear();

      
     

      //cout << "sum rings " <<sumMean4*sumMean4 << " " << sumMean*sumMean<< "  " << sumMean5 +sumMean6<< " sum up down " << sumMean2 << " sum theta "  << sumMean3 << " sum grads " << sumGrads <<endl;

    }
  
    theReturn = theReturn + sumMean*sumMean+sumGrads+sumMean4*sumMean4+sumMean2 + sumMean3 +sumMean6 +sumMean5;

    // theReturn = theReturn + sumMean*sumMean+sumMean4*sumMean4+sumMean6 +sumMean5;

  }

 cout << " " << endl;
       cout << "sum rings " <<sumMean4*sumMean4 << " " << sumMean*sumMean<< "  " << sumMean5 +sumMean6<< " sum up down " << sumMean2 << " sum theta "  << sumMean3 << " sum grads " << sumGrads <<endl;

 canSurf->Update();
   
 delete tempAntGraph;
 delete tempAntGraph2;
 delete tempAntGraph3;
 
 delete myMG;
 delete myMG2;
 delete myMG3;

  return theReturn;
   
}

vector<double> leastSquares(double *xIn, double *yIn, int count) {
  
  double SUMx, SUMy, SUMxy, SUMxx, SUMres, res, slope,
    y_intercept, y_estimate ;
  int i,n;
  

  double x=0;
  double y =0;
  double upper = 0;
  double lower = 100000;
  
  for(int i = 0; i <count; i++){

    if(xIn[i]<lower){
      lower=xIn[i];}
    if(xIn[i]>upper){
      upper=xIn[i];}
  }

  for(int i = 0; i<count; i++){
    if((upper-lower>100)&&xIn[i]<300){
      xIn[i]=xIn[i]+360;
    }  
  }

  SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;
  for (i=0; i<count; i++) {
    
    SUMx = SUMx + xIn[i];
    SUMy = SUMy + yIn[i];
    SUMxy = SUMxy + xIn[i]*yIn[i];
    SUMxx = SUMxx + xIn[i]*xIn[i];
  }
  
  slope = ( SUMx*SUMy - count*SUMxy ) / ( SUMx*SUMx - count*SUMxx );
  y_intercept = ( SUMy - slope*SUMx ) / count;
  
  vector<double> theReturn2;
  theReturn2.push_back(slope);
  theReturn2.push_back(y_intercept);
  return theReturn2;
  
}




