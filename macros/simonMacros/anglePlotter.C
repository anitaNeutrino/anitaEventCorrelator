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


void anglePlotter();

void anglePlotter(){

TCanvas *canSurf= new TCanvas("canSurf","canSurf");
canSurf->Divide(1,3);

   double  deltaTArrayMod[40]={0};
  Double_t deltaR[40]={0};  
   Double_t deltaZ[40]={0};  
   Double_t deltaPhi[40]={0};  
   Double_t deltaHeading[1]={0};

   for(int i = 0;i<40;i++){
     deltaPhi[i] = 0.0;
   }
 
  //TCanvas *canSurf2= new TCanvas("canSurf2","canSurf2");
  // canSurf->Divide(2,3);
 
  double theReturn = 0;
  double sumMean = 0;
  double sumMean2 = 0;
  int count8 = 0;


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
  double dummyArray[40][1] ={{0}}; 
  TGraph *tempAntGraph;

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

  for(int i =0; i < 40; i++){
    phiAngleArray2.push_back(temp);
    deltaTArray2.push_back(temp);

  }

  for(int loop = 1; loop <4; loop++){
    int run = 16+loop;
    //int run = 18;

    canSurf->cd(loop+1);

    phiAngle.clear();
    deltaTVec.clear();
    firstAntVec.clear();
    secondAntVec.clear();

    for(int i = 0; i < 40; i++){
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
     
      pat->pitch=0.64;
      pat->roll=0.14;
      
      pitch=pat->pitch;
      roll=pat->roll;

      UsefulAdu5Pat usefulPat(pat);

      

      for(corInd=0;corInd<28;corInd++) {


	firstAnt=corSum->firstAnt[corInd];
	secondAnt=corSum->secondAnt[corInd];

	if(firstAnt<16){
	  //	  UsefulAdu5Pat usefulPat(pat,-0.042685,0,0);
	  //UsefulAdu5Pat usefulPat(pat,0.0,0,0);
	}else if(firstAnt < 32){
	  //  UsefulAdu5Pat usefulPat(pat,-0.00175,0,0);
	  // if(firstAnt%2){
	  //	      UsefulAdu5Pat usefulPat(pat,0.00653,0,0);
	      // }else{
	      //   UsefulAdu5Pat usefulPat(pat,-0.00175,0,0);
	    //}
	}else{
	  
	  //  UsefulAdu5Pat usefulPat(pat,0.1927,0,0);
	}

	//replace taylor dome
	//	 deltaTExpected=usefulPat.getDeltaTWillySeavey(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	usefulPat.fSourceLongitude=0;
	//	deltaTExpected=usefulPat.getDeltaTTaylor(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
			deltaTExpected=usefulPat.getDeltaTTaylorOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd],deltaR,deltaZ,deltaPhi);

	//	 histSimpleDtDiff->Fill(corSum->maxCorTimes[corInd]-deltaTExpected);
	int antInd=corSum->firstAnt[corInd]%16;
	//	 histDtDiffSep[antInd][labChip][corInd]->Fill(corSum->maxCorTimes[corInd]-deltaTExpected);


	deltaT=corSum->maxCorTimes[corInd];
	maxAnt=corSum->centreAntenna;
	phiWave=usefulPat.getPhiWave()*TMath::RadToDeg();
	phiMaxAnt=fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna)*TMath::RadToDeg();
	corPeak=corSum->maxCorVals[corInd];
	corRMS=corSum->rmsCorVals[corInd];

	Double_t meanPhiAntPair=fGeomTool->getMeanAntPairPhiRelToAftFore(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	Double_t deltaPhiAntPair=fGeomTool->getPhiDiff(meanPhiAntPair,usefulPat.getPhiWave());

	//cout << deltaT << "  "  << deltaTExpected << endl;
	


	if((deltaT - deltaTExpected)*(deltaT - deltaTExpected) < 1 && (corPeak/corRMS)>8 && TMath::Abs(deltaPhiAntPair)*TMath::RadToDeg()<20){
	  phiAngle[0].push_back(phiWave);
	  deltaTVec[0].push_back(deltaT - deltaTExpected + deltaTArrayMod[firstAnt] - deltaTArrayMod[secondAnt]);
	  firstAntVec[0].push_back(firstAnt);
	  secondAntVec[0].push_back(secondAnt);  
	}



      }

      counter++; 
    
    }

 

    double deltaTArray[40][6000] = {{0}};
    double phiAngleArray[40][6000]= {{0}};

    double deltaTArrayCut[40][6000]= {{0}};
    double phiAngleArrayCut[40][6000]= {{0}};

    int middleAnt; 
    int leftAnt,rightAnt;

    int countArray[40] = {0};


    //fill arrays
    //for(int ants = par[0]; ants < par[0]+1; ants++){
    for(int ants = 0; ants < 40; ants++){

      int count = 0;
      int count2 = 0;
      double sumPhi = 0;
      bool true1 = false;
      bool true2 = false;

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

      if(firstAntTemp < 32){
	if((firstAntTemp == ants) &&  (secondAntTemp == rightTemp)){
	  deltaTArray[ants][count] = deltaTVec[0][events];
	  phiAngleArray[ants][count] = phiAngle[0][events];

	  count++;
	}
      }else{

	rightTemp = firstAntTemp+1;
	if(rightTemp>39){
	  rightTemp = 32;
	}

	if(firstAntTemp == ants){
	  
	}

	if((firstAntTemp == ants) &&  (secondAntTemp == rightTemp)){
	  
	  //if(firstAntTemp!=34){
	  //cout << firstAntTemp << "  " << secondAntTemp << endl;
		//    }

	  //cout << ants << "  "<<firstAntTemp << "  " << secondAntTemp << "  " << rightTemp<<endl;

	  deltaTArray[ants][count] = deltaTVec[0][events];
	  phiAngleArray[ants][count] = phiAngle[0][events];
	  
	count++;
	}
      }


      }
	countArray[ants] = count;
      
    }
 

    //make cuts
    for(int ants = 0; ants < 40; ants++){
      int count = 0;
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

 // fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
      double sumPhi = 0;
 


      
     

      for(int events = 0; events < countArray[ants]; events++){

		
		    phiAngleArrayCut[ants][count] = phiAngleArray[ants][events];
		    deltaTArrayCut[ants][count] = deltaTArray[ants][events];
		    
		    count++;
		

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
  for(int ants = 0; ants < 40; ants++){
    count8 = 0;
    for(int events = 1; events < phiAngleArray2[ants].size(); events++){

      if( deltaTArrayLoop[count8]<1){      
      deltaTArrayLoop[count8] = deltaTArray2[ants][events];
      phiAngleArrayLoop[count8] = phiAngleArray2[ants][events];
      count8++;
      }

    }


    if(ants>32){
      //cout << count8 << endl;
    }

    if(count8==0){
      	tempAntGraph  = new TGraph(1, dummyArray[ants], dummyArray[ants]);
    }else{
    
      tempAntGraph  = new TGraph(count8-1,  phiAngleArrayLoop, deltaTArrayLoop);
    

      if(ants <16){  
      
	canSurf->cd(1);

	tempAntGraph->SetMinimum(-0.5);
	tempAntGraph->SetMaximum(0.5);

	tempAntGraph->Draw("ap");    

	tempAntGraph->SetMarkerStyle(1);
   
	tempAntGraph->GetXaxis()->SetLimits(0,360);
    
	


		cout << ants << "  " << count8 <<"  " << tempAntGraph->GetMean(1) << "  " << tempAntGraph->GetMean(2) << endl;

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
    
tempAntGraph->GetXaxis()->SetTitle("phi (degrees)");
tempAntGraph->GetYaxis()->SetTitle("actual - expected time");
	myMG2->Add(tempAntGraph);
	myMG2->Draw("p");


      }else if(ants < 32){
      

	canSurf->cd(2);

	//  canSurf->cd(ants+1);
      
	tempAntGraph->SetMinimum(-0.5);
	tempAntGraph->SetMaximum(0.5);
	tempAntGraph->Draw("ap");
	
	tempAntGraph->SetMarkerStyle(1);
   
	tempAntGraph->GetXaxis()->SetLimits(0,360);
	tempAntGraph->GetXaxis()->SetTitle("phi (degrees)");
    
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
    
tempAntGraph->GetYaxis()->SetTitle("actual - expected time");
	myMG->Add(tempAntGraph);
	myMG->Draw("p");  

      }else {
      
     

	canSurf->cd(3);

	//  canSurf->cd(ants+1);
      
	tempAntGraph->SetMinimum(-0.5);
	tempAntGraph->SetMaximum(0.5);

	double x=0;
	double y =0;
	double upper = 0;
	double lower = 100000;
	for(int i = 0; i <tempAntGraph->GetN(); i++){

	  tempAntGraph->GetPoint(i,x,y);
	  if(x<lower){
	    lower=x;}
	  if(x>upper){
	    upper=x;}
	}





	  for(int i = 0; i <tempAntGraph->GetN(); i++){
	    
	    tempAntGraph->GetPoint(i,x,y);

	    if((upper-lower>100)&&x<300){
	   
	      tempAntGraph->SetPoint(i,x+360,y);
	      lower = 300;
	      upper = 400;   

	    }
	    
	  }


	tempAntGraph->Draw("ap");
	

	
	cout <<"lower " << lower << " " << upper << endl;       

	//	TF1 *fa3 = new TF1("fa3","pol1",lower,upper);
	//tempAntGraph->Fit("fa3","R");

	tempAntGraph->SetMarkerStyle(1);
   
	tempAntGraph->GetXaxis()->SetLimits(0,360);
    
tempAntGraph->GetXaxis()->SetTitle("phi (degrees)");
 tempAntGraph->GetYaxis()->SetTitle("actual - expected time");

 //	 cout << ants <<"  " <<tempAntGraph->GetMean(2) << endl;
	//	sumMean2 = sumMean2 + tempAntGraph->GetMean(2);

	if(ants == 8 || ants == 16 || ants == 12 || ants == 24 || ants == 32){
	  tempAntGraph->SetMarkerColor(8);
	}

	if(ants == 3 || ants == 7 || ants == 23 || ants == 31 || ants == 33){
	  tempAntGraph->SetMarkerColor(1);
	}
    
	if(ants == 9 || ants == 13 || ants == 18 || ants == 26 || ants == 34){
	  tempAntGraph->SetMarkerColor(2);
	}
    
    
	if(ants == 10 || ants == 14 || ants == 20 || ants == 28 || ants == 35){
	  tempAntGraph->SetMarkerColor(3);
	}
    
	if(ants == 11 || ants == 15 || ants == 22 || ants == 30 || ants == 36){
	  tempAntGraph->SetMarkerColor(4);
	}
    
	if(ants == 0 || ants == 4 || ants == 17 || ants == 25 || ants == 37){
	  tempAntGraph->SetMarkerColor(5);
	}


	if(ants == 2 || ants == 6 || ants == 21 || ants == 29 || ants == 38){
	  tempAntGraph->SetMarkerColor(6);
	}

	if(ants == 1 || ants == 5 || ants == 19 || ants == 27 || ants == 39){
	  tempAntGraph->SetMarkerColor(7);
	}
    
	myMG3->Add(tempAntGraph);
	myMG3->Draw("p");  

      }
    } 


  }
 

  cout << sumMean << "  " << sumMean2 <<endl;


   
}

  


