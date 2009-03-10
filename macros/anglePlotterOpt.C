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

double deltaTArrayMod[32];
TCanvas *canSurf= new TCanvas("canSurf","canSurf");
  TMultiGraph *myMG = new TMultiGraph;
  TMultiGraph *myMG3 = new TMultiGraph;
  TMultiGraph *myMG2 = new TMultiGraph;;


void setArray();
void setArray(){
  for(int i = 0; i < 32; i++){
    deltaTArrayMod[i] = 0;
  }

  // deltaTArrayMod[9] = -0.250658;
  // deltaTArrayMod[1] = -0.27446;
  // deltaTArrayMod[10] = -0.10332;
  // deltaTArrayMod[2] = -0.320768;
  // deltaTArrayMod[11] = -0.334935;
  // deltaTArrayMod[3] = -0.197067;
  // deltaTArrayMod[12] = -0.00332853;
  // deltaTArrayMod[4] = -0.194346;
  // deltaTArrayMod[13] = -0.207093;
  // deltaTArrayMod[5] = -0.0954051;
  // deltaTArrayMod[14] = -0.00457999;
  // deltaTArrayMod[6] = -0.144045;
  // deltaTArrayMod[15] = -0.177085;
  // deltaTArrayMod[7] = -0.153101;
  // deltaTArrayMod[8] = 0.0158212;
  // deltaTArrayMod[0] = -0.210793;

}


void setValue(int x, double t);
void setValue(int x, double t){
  
  deltaTArrayMod[x] = t;
  
}

void printArray();
void printArray(){
  
  for(int i = 0; i < 16; i++){
    cout <<  deltaTArrayMod[i] << endl;;
  }
  cout << " " << endl;
}


 
 



double anglePlotterOpt(double *par);  
double anglePlotterOpt(double *par){
  
 



  // double  deltaTArrayMod[32]={0};

  //  TCanvas *canSurf= new TCanvas("canSurf","canSurf");
  //TCanvas *canSurf2= new TCanvas("canSurf2","canSurf2");
  // canSurf->Divide(2,3);
  canSurf->Divide(1,2);
  double theReturn = 0;
  double sumMean = 0;
  double sumMean2 = 0;
  int count8 = 0;
TGraph *tempAntGraph;


  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corrName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  char baseDir[FILENAME_MAX];
  char *corTreeDir = "../../Outfiles";
  double dummyArray[32][1] ={{0}}; 

  vector<vector<double> > phiAngle;
  vector<vector<double> > deltaTVec;
  vector<vector<int> > firstAntVec;
  vector<vector<int> > secondAntVec;

  vector<vector<double> > phiAngleArray2;
  vector<vector<double> > deltaTArray2;

  vector<double> temp;
  vector<int> temp2;
  temp.push_back(0);
  temp2.push_back(0);

  double deltaTArrayLoop[6000] ={0};
  double phiAngleArrayLoop[6000] = {0};


  int leftOpt, rightOpt;
fGeomTool->getThetaPartners(par[0],leftOpt,rightOpt); 
deltaTArrayMod[rightOpt] = par[1]; 
// for(int i = 0; i < 16; i++){
//     cout <<  deltaTArrayMod[i] << endl;;
//   }



   double meanPhi[32] = {0}; 
      meanPhi[0] =22.5;
      meanPhi[1] =67.5;
      meanPhi[2] =112.5;
      meanPhi[3] =157.5;
      meanPhi[4] =202.5;
      meanPhi[5] =247.5;
      meanPhi[6] =292.5;
      meanPhi[7] =337.5;
      meanPhi[8] =45;
      meanPhi[9] =90;
      meanPhi[10] =135;
      meanPhi[11] =180;
      meanPhi[12] =225;
      meanPhi[13] =270;
      meanPhi[14] =315;
      meanPhi[15] =360;

      meanPhi[16] = 22.5;
	meanPhi[17] = 45;
	meanPhi[18] = 67.5;
	meanPhi[19] = 90;
	meanPhi[20] = 112.5;
	meanPhi[21] = 135;
	meanPhi[22] = 157.5;
	meanPhi[23] = 180;
	meanPhi[24] = 202.5;
	meanPhi[25] = 225;
	meanPhi[26] = 247.5;
	meanPhi[27] = 270;
	meanPhi[28] = 292.5;
	meanPhi[29] = 315;
	meanPhi[30] = 337.5;
	meanPhi[31] = 360;




  for(int i = 0; i < 32; i++){
    phiAngleArray2.push_back(temp);
    deltaTArray2.push_back(temp);

  }

  for(int loop = 0; loop < 4; loop++){
    int run = 16+loop;
    //int run = 18;

    canSurf->cd(loop+1);

    phiAngle.clear();
    deltaTVec.clear();
    firstAntVec.clear();
    secondAntVec.clear();

    for(int i = 0; i < 32; i++){
      phiAngle.push_back(temp);
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
      pitch=pat->pitch;
      roll=pat->roll;
      UsefulAdu5Pat usefulPat(pat);
     


      for(corInd=0;corInd<19;corInd++) {


	//replace taylor dome
	//	 deltaTExpected=usefulPat.getDeltaTWillySeavey(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	usefulPat.fSourceLongitude=0;
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
	
	if((deltaT - deltaTExpected)*(deltaT - deltaTExpected) < 1 && (corPeak/corRMS)>8 ){
	  phiAngle[0].push_back(phiWave);
	  deltaTVec[0].push_back(deltaT - deltaTExpected + deltaTArrayMod[firstAnt] - deltaTArrayMod[secondAnt]);
	  firstAntVec[0].push_back(firstAnt);
	  secondAntVec[0].push_back(secondAnt);
	}
      }

      counter++; 
    
    }

 

    double deltaTArray[32][5000] = {{0}};
    double phiAngleArray[32][5000]= {{0}};

    double deltaTArrayCut[32][5000]= {{0}};
    double phiAngleArrayCut[32][5000]= {{0}};

    int middleAnt; 
    int leftAnt,rightAnt;

     int countArray[32];


    //fill arrays
    //for(int ants = par[0]; ants < par[0]+1; ants++){
    for(int ants = 0; ants < 32; ants++){

      int count = 0;
      int count2 = 0;
      double sumPhi = 0;
      bool true1 = false;
      bool true2 = false;

      fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
    
      for(int events = 1; events < phiAngle[0].size(); events++){
	int firstAntTemp = (int)firstAntVec[0][events];
	int secondAntTemp = (int)secondAntVec[0][events];
	int rightTemp = int(rightAnt);

	if((firstAntTemp == ants) &&  (secondAntTemp == rightTemp)){
	  deltaTArray[ants][count] = deltaTVec[0][events];
	  phiAngleArray[ants][count] = phiAngle[0][events];

	  count++;
	}



      }
	countArray[ants] = count;
      
    }

 
    //make cuts
    for(int ants = 0; ants < 32; ants++){
      int count = 0;
      fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
      double sumPhi = 0;
 

      	double lower  = meanPhi[ants]-30;
      	double upper = meanPhi[ants]+30;
      //	double lower  = 0;
      //      	double upper = 360;


	if(lower < 0){
	  lower = 0;
	  upper = 60;
	}

	if(upper > 360){
	  lower = 300;
	  upper = 360;
	}

      //cout <<  ants << "  " << meanPhi[ants] << endl;
	//   cout  << "  " << upper << "  " << lower << endl;

      for(int events = 0; events < countArray[ants]; events++){

		if((phiAngleArray[ants][events] > lower ) && (phiAngleArray[ants][events]< upper)){
      
		 
		    phiAngleArrayCut[ants][count] = phiAngleArray[ants][events];
		    deltaTArrayCut[ants][count] = deltaTArray[ants][events];
		    
		    count++;
		 

	  	}

      }

      for(int events = 0; events < count-1; events++){
	phiAngleArray2[ants].push_back(phiAngleArrayCut[ants][events]);
	deltaTArray2[ants].push_back(deltaTArrayCut[ants][events]);
      }

    }

    delete event;
    delete hk; 
    delete header;
    delete pat;
    delete corSum;

    delete fpHead;
    delete fpGps ;
    delete fpCor;

  }
 

  


  sumMean = 0;
  sumMean2 = 0;
  for(int ants = 0; ants < 32; ants++){
    count8 = 0;
    for(int events = 1; events < phiAngleArray2[ants].size(); events++){

      if( deltaTArrayLoop[count8]<1){      
      deltaTArrayLoop[count8] = deltaTArray2[ants][events];
      phiAngleArrayLoop[count8] = phiAngleArray2[ants][events];
      count8++;
      }

    }

    if(count8==1){
      	tempAntGraph  = new TGraph(1, dummyArray[ants], dummyArray[ants]);
    }else{
    
      tempAntGraph  = new TGraph(count8-1,  phiAngleArrayLoop, deltaTArrayLoop);
    
      if(ants==par[0]){
	theReturn = tempAntGraph->GetMean(2)*tempAntGraph->GetMean(2);
      }


      if(ants <16){  
      
	canSurf->cd(1);

	tempAntGraph->SetMinimum(-1.5);
	tempAntGraph->SetMaximum(1.5);

	tempAntGraph->Draw("ap");    

	tempAntGraph->SetMarkerStyle(1);
   
	tempAntGraph->GetXaxis()->SetLimits(0,360);
    
	//	cout << ants << "  " << count8 <<"  " << tempAntGraph->GetMean(1) << "  " << tempAntGraph->GetMean(2) << endl;

	//	cout << "deltaTArrayMod["<<ants<<"] = " <<  tempAntGraph->GetMean(2) << ";" << endl;

	sumMean = sumMean + tempAntGraph->GetMean(2);

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
    
	//	myMG2->Add(tempAntGraph);
	//	myMG2->Draw("p");
    canSurf->Update();

      }else{
      

	canSurf->cd(2);

	//  canSurf->cd(ants+1);
      
	tempAntGraph->SetMinimum(-1.5);
	tempAntGraph->SetMaximum(1.5);
	tempAntGraph->Draw("ap");
	
	tempAntGraph->SetMarkerStyle(1);
   
	tempAntGraph->GetXaxis()->SetLimits(0,360);
    
	//	 cout << ants <<"  " <<tempAntGraph->GetMean(2) << endl;
	sumMean2 = sumMean2 + tempAntGraph->GetMean(2);

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
    
	//	myMG->Add(tempAntGraph);
	
  canSurf->Update();


      }
    } 
    // delete tempAntGraph;

  }


  //myMG->Draw("ap");  
  //myMG2->Draw("ap");  
 canSurf->Update();
 delete tempAntGraph;

 // cout << sumMean << "  " << sumMean2 <<  "  " <<theReturn <<endl;

  return theReturn;

}

