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




double seavyPlot2();
vector<double> leastSquares(double *xIn, double *yIn, int count);


double seavyPlot2(){

  bool removeTimeOffset = false;
  int corCut = 5; 


  TH1F *calculateMean[40];

  char title[400];
  for(int i = 0; i < 40; i++){
    sprintf(title,"calculateMean%d",i);
    calculateMean[i] = new TH1F(title,title, 100,-1,1); 
  }  

TF1 *f1 = new TF1("f1", "gaus", -1, 1);

  UInt_t eventNumber, triggerTime, triggerTimeNs;
  Int_t firstAnt,secondAnt,maxAnt,corInd;
  Double_t deltaT,deltaTExpected;
  Double_t phiWave, phiMaxAnt;
  Double_t corPeak, corRMS;
  Double_t balloonLat, balloonLon, balloonAlt;
  Double_t heading,pitch,roll;
  Double_t deltaTVec2;    
  Double_t thetaWave;



  int upperAntNums[NUM_PHI]={8,0,9,1,10,2,11,3,12,4,13,5,14,6,15,7};
  int lowerAntNums[NUM_PHI]={16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
  int lowerAntFromUpper[NUM_PHI]={17,19,21,23,25,27,29,31,16,18,20,22,24,26,28,30};
  int nadirAntNums[NUM_PHI]={32,-1,33,-1,34,-1,35,-1,36,-1,37,-1,38,-1,39,-1};
  int bottomFromNadir[8]={16,18,20,22,24,26,28,30};


  Int_t firstAnt2[40]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
  Int_t secondAnt2[40]={9,10,11,12,13,14,15,8,0,1,2,3,4,5,6,7,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,16,33,34,35,36,37,38,39,32};
  Int_t order[40]={0,9,1,10,2,11,3,12,4,13,5,14,6,15,7,8,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};


  Double_t deltaR[40]={0};  
  Double_t deltaZ[40]={0};  
  Double_t deltaPhi[40]={0};  
  Double_t deltaHeading[1]={0};
  double  deltaTArrayMod[40]={0};
  double  deltaTArrayMod2[40]={0};
 
  int count3 = 0;
  int count17 = 0;

 
  int middleAnt; 
  int leftAnt,rightAnt;

  int countArray = 0;
  int countArrayUpper = 0;

  //double penaltyT = 0;
  double penaltyR = 0;
  double penaltyPhi = 0;
  double penaltyZ = 0;

  //double sigmaT = 0.02;
  double sigmaR = 0.05;
  double sigmaZ = 0.1;
  double sigmaPhi = 0.006;


  //   deltaTArrayMod[0] = 0.0966272; deltaR[0] = -0.057775; deltaPhi[0] = -0.00219158; deltaZ[0] = 0.000263412;
//   deltaTArrayMod[0] = -0.128575; deltaR[0] = 0; deltaPhi[ 0] = 0; deltaZ[0] = 0; 
// deltaTArrayMod[1] = -0.0712834; deltaR[1] = 0; deltaPhi[ 1] = 0; deltaZ[1] = 0; 
// deltaTArrayMod[2] = -0.10947; deltaR[2] = 0; deltaPhi[ 2] = 0; deltaZ[2] = 0; 
// deltaTArrayMod[3] = -0.00673809; deltaR[3] = 0; deltaPhi[ 3] = 0; deltaZ[3] = 0; 
// deltaTArrayMod[4] = -0.0586392; deltaR[4] = 0; deltaPhi[ 4] = 0; deltaZ[4] = 0; 
// deltaTArrayMod[5] = -0.107886; deltaR[5] = 0; deltaPhi[ 5] = 0; deltaZ[5] = 0; 
// deltaTArrayMod[6] = -0.109083; deltaR[6] = 0; deltaPhi[ 6] = 0; deltaZ[6] = 0; 
// deltaTArrayMod[7] = -0.125697; deltaR[7] = 0; deltaPhi[ 7] = 0; deltaZ[7] = 0; 
// deltaTArrayMod[8] = -0.0295459; deltaR[8] = 0; deltaPhi[ 8] = 0; deltaZ[8] = 0; 
// deltaTArrayMod[9] = -0.0777603; deltaR[9] = 0; deltaPhi[ 9] = 0.000158312; deltaZ[9] = 0; 
// deltaTArrayMod[10] = -0.171187; deltaR[10] = 0; deltaPhi[ 10] = 0; deltaZ[10] = 0; 
// deltaTArrayMod[11] = -0.138039; deltaR[11] = 0; deltaPhi[ 11] = 0; deltaZ[11] = 0; 
// deltaTArrayMod[12] = -0.0514104; deltaR[12] = 0; deltaPhi[ 12] = 0; deltaZ[12] = 0; 
// deltaTArrayMod[13] = -0.329179; deltaR[13] = 0; deltaPhi[ 13] = 0; deltaZ[13] = 0; 
// deltaTArrayMod[14] = -0.253613; deltaR[14] = 0; deltaPhi[ 14] = 0; deltaZ[14] = 0; 
// deltaTArrayMod[15] = 0.0234044; deltaR[15] = 0; deltaPhi[ 15] = 0; deltaZ[15] = 0; 
// deltaTArrayMod[16] = -0.0195308; deltaR[16] = 0; deltaPhi[ 16] = 0; deltaZ[16] = 0; 
// deltaTArrayMod[17] = -0.0811531; deltaR[17] = 0; deltaPhi[ 17] = 0; deltaZ[17] = 0; 
// deltaTArrayMod[18] = -0.128951; deltaR[18] = 0; deltaPhi[ 18] = 0; deltaZ[18] = 0; 
// deltaTArrayMod[19] = -0.149976; deltaR[19] = 0; deltaPhi[ 19] = 0; deltaZ[19] = 0; 
// deltaTArrayMod[20] = -0.0725489; deltaR[20] = -0.0479999; deltaPhi[ 20] = 0; deltaZ[20] = 0; 
// deltaTArrayMod[21] = -0.105717; deltaR[21] = 0; deltaPhi[ 21] = 0; deltaZ[21] = 0; 
// deltaTArrayMod[22] = -0.136003; deltaR[22] = 0; deltaPhi[ 22] = 0; deltaZ[22] = 0; 
// deltaTArrayMod[23] = -0.123998; deltaR[23] = 0; deltaPhi[ 23] = 0; deltaZ[23] = 0; 
// deltaTArrayMod[24] = -0.0313151; deltaR[24] = -0.000824448; deltaPhi[ 24] = 0; deltaZ[24] = 0; 
// deltaTArrayMod[25] = -0.118466; deltaR[25] = 0; deltaPhi[ 25] = 0; deltaZ[25] = 0; 
// deltaTArrayMod[26] = -0.189091; deltaR[26] = 0; deltaPhi[ 26] = 0; deltaZ[26] = 0; 
// deltaTArrayMod[27] = -0.0718261; deltaR[27] = 0; deltaPhi[ 27] = 0; deltaZ[27] = 0; 
// deltaTArrayMod[28] = -0.280509; deltaR[28] = 0; deltaPhi[ 28] = 0; deltaZ[28] = 0; 
// deltaTArrayMod[29] = -0.17044; deltaR[29] = 0; deltaPhi[ 29] = 0; deltaZ[29] = 0; 
// deltaTArrayMod[30] = -0.101841; deltaR[30] = 0; deltaPhi[ 30] = 0; deltaZ[30] = 0; 
// deltaTArrayMod[31] = -0.126857; deltaR[31] = 0; deltaPhi[ 31] = 0; deltaZ[31] = 0; 
// deltaTArrayMod[32] = -1.13101; deltaR[32] = 0; deltaPhi[ 32] = 0; deltaZ[32] = 0; 
// deltaTArrayMod[33] = -0.938858; deltaR[33] = 0; deltaPhi[ 33] = 0; deltaZ[33] = 0; 
// deltaTArrayMod[34] = -0.981519; deltaR[34] = 0; deltaPhi[ 34] = 0; deltaZ[34] = 0; 
// deltaTArrayMod[35] = -1.0533; deltaR[35] = 0; deltaPhi[ 35] = 0; deltaZ[35] = 0; 
// deltaTArrayMod[36] = -0.96071; deltaR[36] = 0; deltaPhi[ 36] = 0; deltaZ[36] = 0; 
// deltaTArrayMod[37] = -1.00874; deltaR[37] = 0; deltaPhi[ 37] = 0; deltaZ[37] = 0; 
// deltaTArrayMod[38] = -1.01037; deltaR[38] = 0; deltaPhi[ 38] = 0; deltaZ[38] = 0; 
// deltaTArrayMod[39] = -1.02233; deltaR[39] = 0; deltaPhi[ 39] = 0; deltaZ[39] = 0;

//  deltaTArrayMod2[0] = -0.0204477; 
//  deltaTArrayMod2[1] = -0.0197648; 
//  deltaTArrayMod2[2] = -0.0146442; 
//  deltaTArrayMod2[3] = -0.00696469;
//  deltaTArrayMod2[4] = -0.00172871;
//  deltaTArrayMod2[5] = 0.00605453; 
//  deltaTArrayMod2[6] = 0.0176024; 
//  deltaTArrayMod2[7] = -0.0271736;
//  deltaTArrayMod2[8] = -0.0205879;
//  deltaTArrayMod2[9] = -0.0186188;
//  deltaTArrayMod2[10] = -0.013868;
//  deltaTArrayMod2[11] = -0.00849595; 
//  deltaTArrayMod2[12] = -0.00505742; 
//  deltaTArrayMod2[13] = 0.00247081; 
//  deltaTArrayMod2[14] = 0.00989389; 
//  deltaTArrayMod2[15] = -0.048264; 
//  deltaTArrayMod2[16] = -0.0098481;
//  deltaTArrayMod2[17] = -0.00994761;
//  deltaTArrayMod2[18] = -0.0113963; 
//  deltaTArrayMod2[19] = -0.011422; 
//  deltaTArrayMod2[20] = -0.0112273;
//  deltaTArrayMod2[21] = -0.0145946;
//  deltaTArrayMod2[22] = -0.0147203;
//  deltaTArrayMod2[23] = -0.0131724;
//  deltaTArrayMod2[24] = -0.0146237;
//  deltaTArrayMod2[25] = -0.00890171; 
//  deltaTArrayMod2[26] = -0.00593528; 
//  deltaTArrayMod2[27] = -0.00139135; 
//  deltaTArrayMod2[28] = -0.00491395; 
//  deltaTArrayMod2[29] = -0.000256727;
//  deltaTArrayMod2[30] = -0.000421407;
//  deltaTArrayMod2[31] = -0.00799234;
//  deltaTArrayMod2[32] = -0.0141384; 
//  deltaTArrayMod2[33] = -0.0171962; 
//  deltaTArrayMod2[34] = -0.0169355; 
//  deltaTArrayMod2[35] = -0.0156447; 
//  deltaTArrayMod2[36] = -0.0228346; 
//  deltaTArrayMod2[37] = -0.0213988; 
//  deltaTArrayMod2[38] = -0.0284095; 
//  deltaTArrayMod2[39] = -0.0247334; 


  for(int i = 0; i<40;i++){
    //deltaR[i]= deltaR[i] + par[i];   
    //deltaZ[i]=deltaZ[i] + par[i+40];
    //deltaPhi[i]=deltaPhi[i] + par[i+40];
    deltaTArrayMod[i]=deltaTArrayMod2[i] + deltaTArrayMod[i];
   
    //deltaR[i]= deltaR[i] + par[i];   
    //deltaZ[i]= deltaZ[i] + par[i+32];
    //deltaPhi[i]= deltaPhi[i] + par[i+64];
    //    deltaTArrayMod[i] = deltaTArrayMod[i] + par[i+96];

    penaltyR = penaltyR + sqrt((deltaR[i]/sigmaR)*(deltaR[i]/sigmaR));
    penaltyZ = penaltyZ + sqrt((deltaZ[i]/sigmaZ) * (deltaZ[i]/sigmaZ));
    penaltyPhi = penaltyPhi + sqrt((deltaPhi[i]/sigmaPhi)*(deltaPhi[i]/sigmaPhi));
    // penaltyT = penaltyT + sqrt((deltaTArrayMod[i]/sigmaT)*(deltaTArrayMod[i]/sigmaT));

  }

  // for(int i = 0; i<40;i++){
  //   cout << " deltaTArrayMod[" <<i<<"] = " << deltaTArrayMod[i] <<"; deltaR[" << i << "] = "<<deltaR[i] << "; deltaPhi[ " << i << "] = " <<deltaPhi[i] << "; deltaZ[" << i << "] = " << deltaZ[i] << ";"<<endl;
  //  }


  penaltyR = penaltyR/32;
  penaltyZ = penaltyZ/32;
  penaltyPhi = penaltyPhi/32;
  //penaltyT = penaltyT/32;


  TMultiGraph *myMG = new TMultiGraph;
  TMultiGraph *myMG3 = new TMultiGraph;
  TMultiGraph *myMG2 = new TMultiGraph;;

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  fGeomTool->useKurtAnitaIINumbers(1);

  
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char hkName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char corrName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  char baseDir[FILENAME_MAX];
  //char *corTreeDir = "/unix/anita1/rjn/corTrees/justTaylor070509";
  //char *corTreeDir = "/unix/anita1/rjn/corTrees/justTaylorNoDt/";
  char *corTreeDir = "/unix/anita1/swb/seaveyCor/";
  double dummyArray[40] ={0}; 
  TGraph *tempAntGraph;
  TGraph *tempAntGraph2;
  TGraph *tempAntGraph3;

 
 

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






  //  int corFile[35] = {13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 159, 160, 161, 162, 163, 164, 165, 166, 167};
  int corFile[4] = { 17,18,163,165};  
  //  int corFile[3] = { 18, 62,159};

  //58,65,60,162,

  TFile *fpOut;
  fpOut= new TFile("myTree.root","RECREATE");
  TTree *corTree2 = new TTree("corTree2","Tree of Correlation Summaries");
  
  for(int loop = 0; loop <1; loop++){
    //int run = corFile[loop];
    int run = 19;

    
    
    corTree2->Branch("thetaWave",&thetaWave,"thetaWave/D");
    corTree2->Branch("phiWave",&phiWave,"phiWave/D");
    corTree2->Branch("firstAnt",&firstAnt,"firstAnt/I");
    corTree2->Branch("secondAnt",&secondAnt,"secondAnt/I");
    corTree2->Branch("deltaTVec2",&deltaTVec2,"deltaTVec2/D");    

    Long64_t entry=0;
    
    //canSurf->cd(loop+1);
    

    //sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
    sprintf(baseDir,"/unix/anita1/flight0809/root/");
    //sprintf(baseDir,"/Users/simonbevan/Desktop/");
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


    
    for(entry=0;entry<numEntries;entry++) {

      corTree->GetEntry(entry);
      Long64_t headEntry=headTree->GetEntryNumberWithIndex(corSum->eventNumber);
      if(headEntry<0) 
	continue;
      headTree->GetEntry(headEntry);
     
      //if( (header->triggerTimeNs>0.5e6) || (header->triggerTimeNs<0.2e6) )  
      //continue; 

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
      // pat->pitch=0.64;
      // pat->roll=0.14;

      pat->pitch=-0.29;
      pat->roll=0.89;

      pitch=pat->pitch;
      roll=pat->roll;

      UsefulAdu5Pat usefulPat(pat);
      
      for(corInd=0;corInd<29;corInd++) {
	
	firstAnt=corSum->firstAnt[corInd];
	secondAnt=corSum->secondAnt[corInd];

	//replace taylor dome
	
	usefulPat.fSourceLongitude=0;
	//	deltaTExpected=usefulPat.getDeltaTTaylor(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	//deltaTExpected=usefulPat.getDeltaTTaylorOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd],deltaR,deltaZ,deltaPhi);
	deltaTExpected=usefulPat.getDeltaTSeaveyOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd],deltaR,deltaZ,deltaPhi);


	deltaT=corSum->maxCorTimes[corInd];
	maxAnt=corSum->centreAntenna;
	phiWave=usefulPat.getPhiWave()*TMath::RadToDeg();
	thetaWave=usefulPat.getThetaWave()*TMath::RadToDeg();
	phiMaxAnt=fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna)*TMath::RadToDeg();
	corPeak=corSum->maxCorVals[corInd];
	corRMS=corSum->rmsCorVals[corInd];

	Double_t meanPhiAntPair=fGeomTool->getMeanAntPairPhiRelToAftFore(corSum->firstAnt[corInd],corSum->secondAnt[corInd]);
	Double_t deltaPhiAntPair=fGeomTool->getPhiDiff(meanPhiAntPair,usefulPat.getPhiWave());

	//if((corPeak/corRMS)>10 && TMath::Abs(deltaPhiAntPair)*TMath::RadToDeg()<20 && corRMS>40){
	if((corPeak/corRMS)>corCut){
	  //if(corRMS>40){
	  deltaTVec2 = deltaT - deltaTExpected + deltaTArrayMod[firstAnt] - deltaTArrayMod[secondAnt];	 
	  // deltaTVec2 = deltaT - deltaTExpected + deltaTArrayMod[firstAnt] ;	 
	  if(TMath::Abs(deltaTVec2) < 1){
	    corTree2->Fill();  
	  }
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
  }
       
 
  double deltaTChange[40] ={0};

  double theReturn = 0;
  double sumMean = 0;
  double sumMean2 = 0;
  double sumMean3 = 0;
  double sumMean4 = 0;
  double sumMean5 = 0;
  double sumMean6 = 0;
  double sumMean7 = 0;
  double sumMean8 = 0;
  double sumRMS = 0;
  double sumRMS2 = 0;
  double sumRMS3 = 0;
  double sumRMS4 = 0;
  double sumRMS5 = 0;
  int count8 = 0;
  int count18 = 0;
  double sumGrads = 0;

  double deltaTCorrect[40] = {0};

  for(int ants = 0; ants < 40; ants++){
    
    double sumNumbers = 0;
    double meanCount = 0;    

    double sumNumbers2 = 0;
    double meanCount2 = 0;    

    double sumNumbers3 = 0;
    double meanCount3 = 0;    

      
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


    
    for(int events = 0; events < corTree2->GetEntries(); events++){
      
      corTree2->GetEntry(events);
      
      int firstAntTemp = firstAnt;
      int secondAntTemp = secondAnt;
      int rightTemp = int(rightAnt);
      
      int aboveTemp = 0;
      
      if(ants <16){
	aboveTemp = lowerAntFromUpper[ants];
      }else if(ants <32){
	aboveTemp = upperAntNums[ants-16];

      }else{
	aboveTemp = bottomFromNadir[firstAntTemp-32];  
      }
      
      
      if(firstAntTemp < 40){
	
	if( ((firstAntTemp == ants) &&  (secondAntTemp == rightTemp))){
	   
	  if(firstAntTemp<16){
	    sumNumbers = sumNumbers + (double)deltaTVec2;
	    meanCount++;
	  }else if(firstAntTemp<32){
	   
	    sumNumbers2 = sumNumbers2 + (double)deltaTVec2;
	    meanCount2++;
	  }else{
	    sumNumbers3 = sumNumbers3 + (double)deltaTVec2;
	    meanCount3++;
	  }
	  
	  calculateMean[ants]->Fill((float)deltaTVec2);

	}
	


      }
      
      
    }

    
      calculateMean[ants]->Fit("f1", "R");
    if(ants<16){
      if(meanCount == 0){
	deltaTCorrect[ants] = 0;
      }else{
	deltaTCorrect[ants] = f1->GetParameter(1);
	//deltaTCorrect[ants] = sumNumbers/meanCount;
      }
    }else if(ants<32){
      

      if(meanCount2 == 0){
	deltaTCorrect[ants] = 0;
      }else{
	deltaTCorrect[ants] = f1->GetParameter(1);
	//deltaTCorrect[ants] = sumNumbers2/meanCount2;
      }

    }else{
      deltaTCorrect[ants] = f1->GetParameter(1);
      //deltaTCorrect[ants] = sumNumbers3/meanCount3;
    }
 
  }

  double timeOffset[40] = {0};
  for(int j=0;j<40;j++) {
    // timeOffset[j] = deltaTArrayMod[j];
  }


  double sumDeltaTCorrect = 0;
  for(int p = 0; p<16; p++){
    sumDeltaTCorrect = sumDeltaTCorrect + deltaTCorrect[p];
  }
  double meanUpper = sumDeltaTCorrect/16;
  sumDeltaTCorrect = 0;
  
  for(int p = 16; p<32; p++){
    sumDeltaTCorrect = sumDeltaTCorrect + deltaTCorrect[p];
  }
  double meanLower = sumDeltaTCorrect/16;
  sumDeltaTCorrect = 0;  

  for(int p = 32; p<40; p++){
    sumDeltaTCorrect = sumDeltaTCorrect + deltaTCorrect[p];
  }
  double meanNadir = sumDeltaTCorrect/8;
  

  //Here is the main loop
  for(int j=0;j<40;j++) {
    int ant=order[j];
    //Get new diff taking into account whichever antenna t_0 offsets we've already caclulated
    Double_t tempDiff=deltaTCorrect[ant]+(timeOffset[firstAnt2[ant]]-timeOffset[secondAnt2[ant]]);

    //Subtract off the mean of the ring deltaT
    if(ant<16)
      tempDiff-=meanUpper;
    else if(ant<32)
      tempDiff-=meanLower;
    else 
      tempDiff-=meanNadir;
    //And assign the difference to be down to the secondAnt
    timeOffset[secondAnt2[ant]]=tempDiff;
    //    std::cout << secondAnt[ant] << "\t" << antDiff[secondAnt[ant]] << newTimes[ant] << "\n";
  }

  double meanOffSet = TMath::Mean(16,&timeOffset[0]);
  double meanOffSet2 = TMath::Mean(8,&timeOffset[32]);
  
 

  for(int j=0;j<32;j++) {
    timeOffset[j] -= meanOffSet; 
  }
 
  for(int j=32;j<40;j++) {
    timeOffset[j] -= meanOffSet2; 
  }

  if(removeTimeOffset){
    for(int j=0;j<40;j++) {
      timeOffset[j] = 0; 
    }
  }
  double sumNumbers4 = 0;
  double meanCount4 = 0;    

  double sumNumbers5 = 0;
  double meanCount5 = 0;   

  for(int ants = 0; ants < 40; ants++){
    
       
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

    
    for(int events = 0; events < corTree2->GetEntries(); events++){
      
      corTree2->GetEntry(events);
      
      int firstAntTemp = firstAnt;
      int secondAntTemp = secondAnt;
      int rightTemp = int(rightAnt);
      
      int aboveTemp = 0;
      
      if(ants <16){
	aboveTemp = lowerAntFromUpper[ants];
      }else if(ants <32){
	aboveTemp = upperAntNums[ants-16];
      }else{
	aboveTemp = bottomFromNadir[firstAntTemp-32];  
      }
      
      
      if(firstAntTemp < 40){
	
	if(((firstAntTemp == ants) &&  (secondAntTemp == aboveTemp))){
	  

	  if(firstAntTemp<32){
	    sumNumbers4 = sumNumbers4+ (double)deltaTVec2 + timeOffset[firstAntTemp] - timeOffset[secondAntTemp];
	    meanCount4++;
	  }else{
	    sumNumbers5 = sumNumbers5+ (double)deltaTVec2 + timeOffset[firstAntTemp] - timeOffset[secondAntTemp];
	    meanCount5++;
	  } 
	}
      }
    }
  }


  double  deltaTCorrect2 = sumNumbers4/meanCount4;
  if(meanCount4==0){
    deltaTCorrect2 = 0;
  }

  double   deltaTCorrect3 = sumNumbers5/meanCount5;
  if(meanCount5==0){
    deltaTCorrect3 = 0;
  }


  //cout << deltaTCorrect2  << "  "<< deltaTCorrect3 << endl;

  for(int ants = 0; ants < 40; ants++){
    
    if(ants < 16){
      timeOffset[ants] -= deltaTCorrect2; 
    }
    
    if(ants>31){
      timeOffset[ants] -= deltaTCorrect3; 
    }
    
    //   cout << timeOffset[ants] << endl;

    if(removeTimeOffset){
      timeOffset[ants] = 0;
    }

  }




  for(int ants = 0; ants < 40; ants++){

 
    int count = 0;
    int count2 = 0;
    int count16 = 0;
    count17 = 0;
    count3 = 0;
     
    if(ants <40){
      fGeomTool->getThetaPartners(ants,leftAnt,rightAnt); 
    }else{
      leftAnt = ants;
      rightAnt = ants +1;
      if(ants == 39){
	leftAnt = ants;
	rightAnt = 32;
      }
    }    

 
    for(int events = 0; events < corTree2->GetEntries(); events++){

      corTree2->GetEntry(events);

      int firstAntTemp = firstAnt;
      int secondAntTemp = secondAnt;
      int rightTemp = int(rightAnt);

      int aboveTemp = 0;

      if(ants <16){
	aboveTemp = lowerAntFromUpper[ants];
      }else if(ants < 32){
	aboveTemp = upperAntNums[ants-16];
      }else{
	aboveTemp = bottomFromNadir[firstAntTemp-32];  
      }
  
	
      if(firstAntTemp < 40){
      
	if( ((firstAntTemp == ants) &&  (secondAntTemp == rightTemp))){
	 
	  phiAngleArray2.push_back((double)phiWave);
	  

	  if(firstAntTemp < 32){
	  

	    //deltaTArray2.push_back((double)deltaTVec2);
	    deltaTArray2.push_back((double)deltaTVec2 + timeOffset[firstAntTemp] - timeOffset[secondAntTemp]);
	  }else{
	    
	    //  deltaTArray2.push_back((double)deltaTVec2 + deltaTArrayMod[firstAntTemp] - deltaTArrayMod[secondAntTemp]);
	    deltaTArray2.push_back((double)deltaTVec2 + timeOffset[firstAntTemp]  - timeOffset[secondAntTemp]);
	  }


	  //  deltaTArray2.push_back((double)deltaTVec2);
	  
	  
	  // cout << phiWave << endl;

	  count++;

	} else if(((firstAntTemp == ants) &&  (secondAntTemp == aboveTemp))){

	  thetaAngleArray2Upper.push_back((double)thetaWave);
	  phiAngleArray2Upper.push_back((double)phiWave);
	  
	  
	  if(firstAntTemp < 32){
	    //deltaTArray2Upper.push_back((double)deltaTVec2 );
	    deltaTArray2Upper.push_back((double)deltaTVec2 + timeOffset[firstAntTemp] - timeOffset[secondAntTemp]);
	  }else{
	    // deltaTArray2Upper.push_back((double)deltaTVec2 +  deltaTArrayMod[firstAntTemp] -  deltaTArrayMod[secondAntTemp]);
	    deltaTArray2Upper.push_back((double)deltaTVec2 + timeOffset[firstAntTemp]  - timeOffset[secondAntTemp]);
	  }
	
	  count16++;   
	    
	}
 
      }
    
    
    }


    countArray = count;
    countArrayUpper = count16;
  
 
	
    count8 = 0;
    count18 = 0;


    // cout << phiAngleArray2.size() << endl;

   

    double deltaTArrayLoopUpper[10000] = {0};
    double phiAngleArrayLoopUpper[10000] = {0};
    double thetaAngleArrayLoopUpper[10000] = {0};

    double deltaTArrayLoop[10000] ={0};
    double phiAngleArrayLoop[10000] = {0};

    for(int events = 1; events < phiAngleArray2.size(); events++){

      if( deltaTArray2[events]*deltaTArray2[events]  <100){  
	deltaTArrayLoop[count8] = deltaTArray2[events];
	phiAngleArrayLoop[count8] = phiAngleArray2[events];
	count8++;
      
	// cout <<  phiAngleArray2[events] << "  " <<  deltaTArray2[events]<< endl;

      }
      
    }

    // deltaTArray2.clear();
    //  phiAngleArray2.clear();

    for(int events = 1; events < phiAngleArray2Upper.size(); events++){

      if( deltaTArray2Upper[events]* deltaTArray2Upper[events] <100){      
	deltaTArrayLoopUpper[count18] = deltaTArray2Upper[events];
	phiAngleArrayLoopUpper[count18] = phiAngleArray2Upper[events];
	thetaAngleArrayLoopUpper[count18] = thetaAngleArray2Upper[events];
	count18++;
      
      }
      
    }


    //    cout << ants << count << "  " << count16 << endl;
	
    //   cout << count8 << endl;
       

    if(count8==0){

      tempAntGraph  = new TGraph(1, dummyArray, dummyArray);
      tempAntGraph2  = new TGraph(1, dummyArray, dummyArray);
      tempAntGraph3  = new TGraph(1, dummyArray, dummyArray);
    }else{
      tempAntGraph  = new TGraph(count8-1,  phiAngleArrayLoop, deltaTArrayLoop);
      tempAntGraph2  = new TGraph(count18-1,  phiAngleArrayLoopUpper, deltaTArrayLoopUpper);      
      tempAntGraph3  = new TGraph(count18-1,  thetaAngleArrayLoopUpper, deltaTArrayLoopUpper);      

      if(ants < 40){  
	
	//	canSurf->cd(1);
       
	TCanvas *canSurf= new TCanvas("canSurf","canSurf");
	tempAntGraph->SetMinimum(-0.5);
	tempAntGraph->SetMaximum(0.5);
	tempAntGraph->Draw("ap");    

	tempAntGraph->SetMarkerStyle(1);
	tempAntGraph->GetXaxis()->SetLimits(0,360);
    
	

	if(ants < 16){
	  sumMean = sumMean + tempAntGraph->GetMean(2);
	  sumMean5 = sumMean5 + tempAntGraph->GetMean(2)*tempAntGraph->GetMean(2);
	  sumRMS = sumRMS + tempAntGraph->GetRMS(2);
	    
	  //  cout << " leftRight " << count8 << " " << ants << "  " << tempAntGraph->GetMean(2) << "  " << tempAntGraph->GetRMS(2) << endl;

	}else if(ants <32){

	  sumMean4 = sumMean4 + tempAntGraph->GetMean(2);
	  sumMean6 = sumMean6 + tempAntGraph->GetMean(2)*tempAntGraph->GetMean(2);
	  sumRMS2 = sumRMS2 + tempAntGraph->GetRMS(2);

	  // cout <<" leftRight "<< ants << "  " << tempAntGraph->GetMean(2) <<  "  " << tempAntGraph->GetRMS(2)<<endl;

	   
	}else{

	  sumMean7 = sumMean7 + tempAntGraph->GetMean(2);
	  sumMean8 = sumMean8 + tempAntGraph->GetMean(2)*tempAntGraph->GetMean(2);
	  sumRMS5 = sumRMS5 + tempAntGraph->GetRMS(2);

	  // cout << " leftRight "<< ants << "  " << tempAntGraph->GetMean(2) << "  " << tempAntGraph->GetRMS(2) << endl;


	}



	

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


	if(ants == 2 || ants == 6 || ants == 21 || ants == 29 || ants ==38 ){
	  tempAntGraph->SetMarkerColor(6);
	}

	if(ants == 1 || ants == 5 || ants == 19 || ants == 27 || ants == 39){
	  tempAntGraph->SetMarkerColor(7);
	}
    
	
	tempAntGraph->GetXaxis()->SetTitle("phi (degrees)");
	tempAntGraph->GetYaxis()->SetTitle("actual - expected time");
	myMG->Add(tempAntGraph);
	myMG->Draw("p");

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
	TCanvas *canSurf2= new TCanvas("canSurf2","canSurf2");
	tempAntGraph2  = new TGraph(count18-1,  phiAngleArrayLoopUpper, deltaTArrayLoopUpper);
      
	

	tempAntGraph2->SetMinimum(-0.5);
	tempAntGraph2->SetMaximum(0.5);
	tempAntGraph2->Draw("ap");    

	tempAntGraph2->SetMarkerStyle(1);
	tempAntGraph2->GetXaxis()->SetLimits(0,360);
    

	//	cout << " updown " << count18 << "  " <<ants << "  " << tempAntGraph->GetMean(2) << "  " << tempAntGraph->GetRMS(2)<< endl;


	sumMean2 = sumMean2 + tempAntGraph2->GetMean(2)*tempAntGraph2->GetMean(2);
	sumRMS3 = sumRMS3 + tempAntGraph2->GetRMS(2);

	if(ants == 8 || ants == 16 || ants == 12 || ants == 24 || ants == 32){
	  tempAntGraph2->SetMarkerColor(8);
	}

	if(ants == 3 || ants == 7 || ants == 23 || ants == 31 || ants == 33){
	  tempAntGraph2->SetMarkerColor(1);
	}
    
	if(ants == 9 || ants == 13 || ants == 18 || ants == 26 || ants == 34){
	  tempAntGraph2->SetMarkerColor(2);
	}
    
    
	if(ants == 10 || ants == 14 || ants == 20 || ants == 28 || ants == 35){
	  tempAntGraph2->SetMarkerColor(3);
	}
    
	if(ants == 11 || ants == 15 || ants == 22 || ants == 30 || ants == 36){
	  tempAntGraph2->SetMarkerColor(4);
	}
    
	if(ants == 0 || ants == 4 || ants == 17 || ants == 25 || ants == 37){
	  tempAntGraph2->SetMarkerColor(5);
	}


	if(ants == 2 || ants == 6 || ants == 21 || ants == 29 || ants == 38){
	  tempAntGraph2->SetMarkerColor(6);
	}

	if(ants == 1 || ants == 5 || ants == 19 || ants == 27 || ants == 39){
	  tempAntGraph2->SetMarkerColor(7);
	}
    
	tempAntGraph2->GetXaxis()->SetTitle("phi (degrees)");
	tempAntGraph2->GetYaxis()->SetTitle("actual - expected time");
	myMG2->Add(tempAntGraph2);
	myMG2->Draw("p");

	TCanvas *canSurf3= new TCanvas("canSurf3","canSurf3");      
	tempAntGraph3  = new TGraph(count18-1,  thetaAngleArrayLoopUpper, deltaTArrayLoopUpper);
	tempAntGraph3->SetMinimum(-0.5);
	tempAntGraph3->SetMaximum(0.5);
	tempAntGraph3->Draw("ap");    

	tempAntGraph3->SetMarkerStyle(1);
	tempAntGraph3->GetXaxis()->SetLimits(0,360);
    
	sumMean3 = sumMean3 + tempAntGraph3->GetMean(2)*tempAntGraph3->GetMean(2);
	sumRMS4 = sumRMS4 + tempAntGraph3->GetRMS(2);

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
    
	tempAntGraph3->GetXaxis()->SetTitle("theta (degrees)");
	tempAntGraph3->GetYaxis()->SetTitle("actual - expected time");
	myMG3->Add(tempAntGraph3);
	myMG3->Draw("p");


	vector<double> myFit = leastSquares(phiAngleArrayLoopUpper, deltaTArrayLoopUpper, count18-1);
	double slope = myFit[0];
	double intercept = myFit[1];	
       
	sumGrads = sumGrads + myFit[0]*myFit[0];


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


  
  theReturn = theReturn + sumMean*sumMean+sumGrads+sumMean4*sumMean4+ sumMean7*sumMean7+sumMean2 + sumMean3 +sumMean6 +sumMean5 +sumMean8+ sumRMS + sumRMS2 +sumRMS3 + sumRMS4+sumRMS5; //+ penaltyR + penaltyZ + penaltyPhi;

  //theReturn = theReturn + sumMean*sumMean+sumGrads+sumMean4*sumMean4 +sumMean2 + sumMean3 +sumMean6 +sumMean5 + sumRMS + sumRMS2 +sumRMS3 + sumRMS4; //+ penaltyR + penaltyZ + penaltyPhi;

  // theReturn = theReturn + sumMean*sumMean+sumMean4*sumMean4+sumMean6 +sumMean5;

  //corTree2->Print("all");
  //corTree2->AutoSave();
  delete corTree2;   
  delete fpOut;

  // canSurf->Update();

  cout << " " << endl;
  cout << "sum rings " <<sumMean*sumMean << " " << sumMean4*sumMean4<< "  " << sumMean7*sumMean7 << " " <<sumMean5 +sumMean6 + sumMean8<< " sum up down " << sumMean2 << " sum theta "  << sumMean3 << " sum grads " << sumGrads << " RMS sums " << sumRMS<<" " << " " << sumRMS2 << " " << sumRMS5<< "  "<< sumRMS3 << " " << sumRMS4 << endl;
  




  //delete tempAntGraph;
  //delete tempAntGraph2;
  //delete tempAntGraph3;
 

  //delete myMG;
  //delete myMG2;
  //delete myMG3;


  for(int i = 0; i < 16; i++){
    //  timeOffset[i] = timeOffset[i] - deltaTCorrect2;
  }

  for(int i =32; i <40; i++){
    //  timeOffset[i] = timeOffset[i] - deltaTCorrect3;
  }

  for(int i = 0; i<40;i++){
    cout << " deltaTArrayMod[" <<i<<"] = " << deltaTArrayMod[i] <<"; deltaR[" << i << "] = "<<deltaR[i] << "; deltaPhi[ " << i << "] = " <<deltaPhi[i] << "; deltaZ[" << i << "] = " << deltaZ[i] << ";"<<endl;
  }
  cout << deltaTCorrect2  << endl;

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




