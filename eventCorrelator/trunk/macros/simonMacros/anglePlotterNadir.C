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



void anglePlotterNadir();

  
void anglePlotterNadir(){

  double  deltaTArrayMod[40]={0};

deltaTArrayMod[9] = -0.0190765;
deltaTArrayMod[1] = -0.0605068;
deltaTArrayMod[10] = 0.107009;
deltaTArrayMod[2] = -0.122117;
deltaTArrayMod[11] = -0.119334;
deltaTArrayMod[3] = -0.00787629;
deltaTArrayMod[12] = 0.196738;
deltaTArrayMod[4] = -0.0221204;
deltaTArrayMod[13] = -0.040719;
deltaTArrayMod[5] = 0.0509295;
deltaTArrayMod[14] = 0.149295;
deltaTArrayMod[6] = 0.0178677;
deltaTArrayMod[15] = 0.00872976;
deltaTArrayMod[7] = 0.0156516;
deltaTArrayMod[8] = 0.226635;
deltaTArrayMod[0] = -0.000820606;

deltaTArrayMod[17] = -0.0530281;
deltaTArrayMod[18] = -0.14337;
deltaTArrayMod[19] = -0.131643;
deltaTArrayMod[20] = -0.369562;
deltaTArrayMod[21] = -0.323906;
deltaTArrayMod[22] = -0.333173;
deltaTArrayMod[23] = -0.238455;
deltaTArrayMod[24] = -0.167277;
deltaTArrayMod[25] = -0.0935844;
deltaTArrayMod[26] = -0.0914406;
deltaTArrayMod[27] = 0.0895252;
deltaTArrayMod[28] = -0.0376108;
deltaTArrayMod[29] = 0.0255145;
deltaTArrayMod[30] = 0.0604195;
deltaTArrayMod[31] = 0.0288773;
deltaTArrayMod[16] = -0.013401;

 deltaTArrayMod[32] = 0;
 deltaTArrayMod[33] = 0;
 deltaTArrayMod[34] = 0;
 deltaTArrayMod[35] = 0;
 deltaTArrayMod[36] = 0;
 deltaTArrayMod[37] = 0;
 deltaTArrayMod[38] = 0;
 deltaTArrayMod[39] = 0;


 Double_t deltaR[40]={0};  
 Double_t deltaZ[40]={0};  
 Double_t deltaPhi[40]={0};  
 Double_t deltaHeading[1]={0};

 //r
// deltaR[0+32] = -0.164323;
// deltaR[1+32] = -0.0807173;
// deltaR[2+32] = -0.0775714;
// deltaR[3+32] = 0.000727919;
// deltaR[4+32] = 0.00848306;
// deltaR[5+32] = 0.0139346;
// deltaR[6+32] = -0.12191;
// deltaR[7+32] = -0.0196709;
// deltaZ[0+32] = -0.0485368;
// deltaZ[1+32] = 0.0400569;
// deltaZ[2+32] = -0.079289;
// deltaZ[3+32] = -0.00955065;
// deltaZ[4+32] = -0.0416788;
// deltaZ[5+32] = 0.0518347;
// deltaZ[6+32] = -0.0888914;
// deltaZ[7+32] = 0.186359;
// deltaPhi[0+32] = 0;
// deltaPhi[1+32] = 0;
// deltaPhi[2+32] = 0;
// deltaPhi[3+32] = 0;
// deltaPhi[4+32] = 0;
// deltaPhi[5+32] = 0;
// deltaPhi[6+32] = 0;
// deltaPhi[7+32] = 0;

 //r, z
// deltaR[0+32] = -0.0152755;
// deltaR[1+32] = -0.00553843;
// deltaR[2+32] = 0.00212581;
// deltaR[3+32] = 0.00585106;
// deltaR[4+32] = 0.0119651;
// deltaR[5+32] = 0.00752331;
// deltaR[6+32] = -0.00823993;
// deltaR[7+32] = 0.00198857;
// deltaZ[0+32] = 0.0323779;
// deltaZ[1+32] = 0.0136696;
// deltaZ[2+32] = -0.00133407;
// deltaZ[3+32] = -0.0118401;
// deltaZ[4+32] = -0.0263063;
// deltaZ[5+32] = -0.0178185;
// deltaZ[6+32] = 0.0144546;
// deltaZ[7+32] = -0.00283883;
// deltaPhi[0+32] = 0.0024448;
// deltaPhi[1+32] = -0.0220422;
// deltaPhi[2+32] = -0.0238489;
// deltaPhi[3+32] = -0.0133056;
// deltaPhi[4+32] = 0.00765958;
// deltaPhi[5+32] = 0.0325572;
// deltaPhi[6+32] = 0.033349;
// deltaPhi[7+32] = 0.0199525;


//r, z, phi
// deltaR[0+32] = 0.156734;
// deltaR[1+32] = 0.2125;
// deltaR[2+32] = 0.169629;
// deltaR[3+32] = 0.141297;
// deltaR[4+32] = 0.127732;
// deltaR[5+32] = 0.17145;
// deltaR[6+32] = 0.0990739;
// deltaR[7+32] = 0.203309;


// deltaZ[0+32] = 0.163051;
// deltaZ[1+32] = 0.281998;
// deltaZ[2+32] = 0.0430088;
// deltaZ[3+32] = -0.170058;
// deltaZ[4+32] = -0.297955;
// deltaZ[5+32] = -0.114614;
// deltaZ[6+32] = -0.0854255;
// deltaZ[7+32] = 0.14983;


// deltaPhi[0+32] = -0.0107883;
// deltaPhi[1+32] = 0.00241582;
// deltaPhi[2+32] = 0.0254493;
// deltaPhi[3+32] = 0.020575;
// deltaPhi[4+32] = 0.00825393;
// deltaPhi[5+32] = -0.0109487;
// deltaPhi[6+32] = 0.0112327;
// deltaPhi[7+32] = 0.0126042;


//r, phi, deltaT
deltaR[0+32] = 0.192096;
deltaR[1+32] = 0.185281;
deltaR[2+32] = 0.172963;
deltaR[3+32] = 0.136439;
deltaR[4+32] = 0.126339;
deltaR[5+32] = 0.143903;
deltaR[6+32] = 0.144807;
deltaR[7+32] = 0.171424;
deltaTArrayMod[0+32] = -0.33006;
deltaTArrayMod[1+32] = -0.31323;
deltaTArrayMod[2+32] = -0.0845225;
deltaTArrayMod[3+32] = 0.234214;
deltaTArrayMod[4+32] = 0.382551;
deltaTArrayMod[5+32] = 0.232483;
deltaTArrayMod[6+32] = -0.0292506;
deltaTArrayMod[7+32] = -0.122308;
deltaPhi[0+32] = -0.010578;
deltaPhi[1+32] = 0.0020203;
deltaPhi[2+32] = 0.0283709;
deltaPhi[3+32] = 0.0171348;
deltaPhi[4+32] = 0.00764585;
deltaPhi[5+32] = -0.0151912;
deltaPhi[6+32] = 0.00774264;
deltaPhi[7+32] = 0.0128789;


  TCanvas *canSurf= new TCanvas("canSurf","canSurf");
  //TCanvas *canSurf2= new TCanvas("canSurf2","canSurf2");
  // canSurf->Divide(2,3);
  
  double theReturn = 0;
  double sumMean = 0;
  double sumMean2 = 0;
  double sumGrads = 0;
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
  char *corTreeDir = "../../Outfiles";
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

// 	meanPhi[32] = 45-12.5;
// 	meanPhi[33] = 90-12.5;
// 	meanPhi[34] = 135-12.5;
// 	meanPhi[35] = 180-12.5;
// 	meanPhi[36] = 225-12.5;
// 	meanPhi[37] = 270-12.5;
// 	meanPhi[38] = 315-12.5;
// 	meanPhi[39] = 360-12.5;


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
      pitch=pat->pitch;
      roll=pat->roll;
      UsefulAdu5Pat usefulPat(pat);

      

      for(corInd=19;corInd<28;corInd++) {


	firstAnt=corSum->firstAnt[corInd];
	secondAnt=corSum->secondAnt[corInd];

	 UsefulAdu5Pat usefulPat(pat,0.1927,0,0);

	//    UsefulAdu5Pat usefulPat(pat,0.115,0,0);

	    //       UsefulAdu5Pat usefulPat(pat,0,0,0);

	//replace taylor dome
	//	 deltaTExpected=usefulPat.getDeltaTWillySeavey(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	usefulPat.fSourceLongitude=0;
	//deltaTExpected=usefulPat.getDeltaTTaylor(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
 deltaTExpected=usefulPat.getDeltaTTaylorOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd],deltaR,deltaZ,deltaPhi);



	//cout << (deltaT - deltaTExpected + deltaTArrayMod[firstAnt] - deltaTArrayMod[secondAnt]) << endl;
	//	cout << firstAnt << "  "  << secondAnt << "  " <<deltaT<< "  "  << deltaTExpected << "  "  << deltaTArrayMod[firstAnt] << " "  << deltaTArrayMod[secondAnt]<< endl;

	//	 histSimpleDtDiff->Fill(corSum->maxCorTimes[corInd]-deltaTExpected);
	int antInd=corSum->firstAnt[corInd]%16;
	//	 histDtDiffSep[antInd][labChip][corInd]->Fill(corSum->maxCorTimes[corInd]-deltaTExpected);


	deltaT=corSum->maxCorTimes[corInd];
	maxAnt=corSum->centreAntenna;
	phiWave=usefulPat.getPhiWave()*TMath::RadToDeg();
	phiMaxAnt=fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna)*TMath::RadToDeg();
	corPeak=corSum->maxCorVals[corInd];
	corRMS=corSum->rmsCorVals[corInd];
	
		if((deltaT - deltaTExpected)*(deltaT - deltaTExpected) < 100 && (corPeak/corRMS)>8 ){
	//	if((deltaT - deltaTExpected)*(deltaT - deltaTExpected) < 100 ){
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
    for(int antLoop = 0; antLoop < 8; antLoop++){

int      ants = antLoop + 32;

      int count = 0;
      int count2 = 0;
      double sumPhi = 0;
      bool true1 = false;
      bool true2 = false;

      //if(ants <32){
      fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
      //}else{
	//leftAnt = ants;
	//rightAnt = ants +1;
	//if(ants == 39){
      // leftAnt = ants;
      //rightAnt = 32;
      //}
      //}    

   
      for(int events = 1; events < phiAngle[0].size(); events++){
	int firstAntTemp = (int)firstAntVec[0][events];
	int secondAntTemp = (int)secondAntVec[0][events];
	int rightTemp = int(rightAnt);

  

	rightTemp = firstAntTemp+1;
	if(rightTemp>39){
	  rightTemp = 32;
	}



	if((firstAntTemp == ants) &&  (secondAntTemp == rightTemp)){
	 
	  deltaTArray[ants][count] = deltaTVec[0][events];
	  phiAngleArray[ants][count] = phiAngle[0][events];
	  
	count++;
	}
      


      }
	countArray[ants] = count;
      
    }
 

    //make cuts
    for(int antLoop = 0; antLoop < 8; antLoop++){
     int ants = antLoop +32;
      int count = 0;
   
	leftAnt = ants - 1;
	rightAnt = ants + 1;
	if(ants == 39){
	  leftAnt = ants - 1;
	rightAnt = 32;
	
      }    

 // fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
      double sumPhi = 0;
 

           	double lower  = meanPhi[ants]-20;
            	double upper = meanPhi[ants]+20;
		//double lower  = 0;
		//	double upper = 360;


	if(lower < 0){
	  lower = 0;
	  upper = 30;
	}

	if(upper > 360){
	  lower = 330;
	  upper = 360;
	}

      //cout <<  ants << "  " << meanPhi[ants] << endl;
	//   cout  << "  " << upper << "  " << lower << endl;

	// cout << countArray[ants] << endl;

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
  sumGrads = 0;
  for(int antLoop = 0; antLoop < 8; antLoop++){
    count8 = 0;
    int ants = antLoop + 32;

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
    



	//  canSurf->cd(ants+1);
      
	tempAntGraph->SetMinimum(-1.5);
	tempAntGraph->SetMaximum(1.5);
	tempAntGraph->Draw("ap");
	
	tempAntGraph->SetMarkerStyle(1);
   
	tempAntGraph->GetXaxis()->SetLimits(0,360);
    
	//	 cout << ants <<"  " <<tempAntGraph->GetMean(2) << endl;
	//sumMean2 = sumMean2 + tempAntGraph->GetMean(2)* tempAntGraph->GetMean(2);
	sumMean2 = sumMean2 + tempAntGraph->GetMean(2);

		TF1 *f1 = new TF1("f1","[0]*x+[1]",meanPhi[ants]-20,meanPhi[ants]+20);
		tempAntGraph->Fit("f1","R");
		double fit = (f1->GetParameter(0))*100;
		
		sumGrads = sumGrads + fit*fit;
		f1->Draw("same");


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
 
  cout << sumGrads << "  " << sumMean2 <<endl;
   
}

