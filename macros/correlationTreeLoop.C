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
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>



void correlationTreeLoop(int run, char *baseDir, char *corTreeDir, char *outputDir);


double antTimeOffset[40];
double deltaR[40];
double deltaPhi[40];
double deltaZ[40];

void correlationTreeLoop(int run,char *baseDir, char *corTreeDir, char *outputDir) {
// antTimeOffset[0] = 0.292591; deltaR[0] = -0.0384839; deltaPhi[ 0] = -0.0100608; deltaZ[0] = 0;
// antTimeOffset[1] = 0.172375; deltaR[1] = 0.00634697; deltaPhi[ 1] = -0.00313443; deltaZ[1] = 0;
// antTimeOffset[2] = 0.342454; deltaR[2] = -0.0861167; deltaPhi[ 2] = -0.015312; deltaZ[2] = 0;
// antTimeOffset[3] = 0.0248334; deltaR[3] = 0.0461873; deltaPhi[ 3] = 0.00206827; deltaZ[3] = 0;
// antTimeOffset[4] = 0.0372539; deltaR[4] = 0.0153388; deltaPhi[ 4] = -0.0227948; deltaZ[4] = 0;
// antTimeOffset[5] = 0.137506; deltaR[5] = -0.00927728; deltaPhi[ 5] = 0.00750385; deltaZ[5] = 0;
// antTimeOffset[6] = 0.0475841; deltaR[6] = 0.0239867; deltaPhi[ 6] = 0.00388065; deltaZ[6] = 0;
// antTimeOffset[7] = 0.123741; deltaR[7] = 0.0125282; deltaPhi[ 7] = -0.00131021; deltaZ[7] = 0;
// antTimeOffset[8] = 0.380551; deltaR[8] = -0.0111636; deltaPhi[ 8] = -0.0299233; deltaZ[8] = 0;
// antTimeOffset[9] = 0.445956; deltaR[9] = -0.0959452; deltaPhi[ 9] = -0.00165365; deltaZ[9] = 0;
// antTimeOffset[10] = 0.439757; deltaR[10] = -0.0330808; deltaPhi[ 10] = -0.0107407; deltaZ[10] = 0;
// antTimeOffset[11] = 0.244031; deltaR[11] = -0.0475617; deltaPhi[ 11] = 0.0145914; deltaZ[11] = 0;
// antTimeOffset[12] = 0.321419; deltaR[12] = 0.0196292; deltaPhi[ 12] = -0.0150373; deltaZ[12] = 0;
// antTimeOffset[13] = 0.118473; deltaR[13] = -0.0190837; deltaPhi[ 13] = -0.0121967; deltaZ[13] = 0;
// antTimeOffset[14] = 0.27363; deltaR[14] = -0.00922367; deltaPhi[ 14] = -0.0038106; deltaZ[14] = 0;
// antTimeOffset[15] = 0.187377; deltaR[15] = -0.0294811; deltaPhi[ 15] = 0.0106842; deltaZ[15] = 0;
// antTimeOffset[16] = 0.0682454; deltaR[16] = 0.0140245; deltaPhi[ 16] = -0.0087849; deltaZ[16] = 0;
// antTimeOffset[17] = 0.296899; deltaR[17] = -0.0621836; deltaPhi[ 17] = 0.000682206; deltaZ[17] = 0;
// antTimeOffset[18] = 0.180588; deltaR[18] = -0.0379325; deltaPhi[ 18] = -0.00516052; deltaZ[18] = 0;
// antTimeOffset[19] = 0.165257; deltaR[19] = -0.0108062; deltaPhi[ 19] = -0.00770935; deltaZ[19] = 0;
// antTimeOffset[20] = 0.176423; deltaR[20] = -0.0601935; deltaPhi[ 20] = -0.00862535; deltaZ[20] = 0;
// antTimeOffset[21] = 0.301777; deltaR[21] = -0.0968276; deltaPhi[ 21] = -0.00920648; deltaZ[21] = 0;
// antTimeOffset[22] = 0.0969425; deltaR[22] = -0.0348523; deltaPhi[ 22] = 0.00037431; deltaZ[22] = 0;
// antTimeOffset[23] = 0.0422886; deltaR[23] = 0.0121726; deltaPhi[ 23] = 0.00310935; deltaZ[23] = 0;
// antTimeOffset[24] = -0.01737; deltaR[24] = 0.0405193; deltaPhi[ 24] = -0.00546085; deltaZ[24] = 0;
// antTimeOffset[25] = 0.0412981; deltaR[25] = 0.0239992; deltaPhi[ 25] = -0.00901249; deltaZ[25] = 0;
// antTimeOffset[26] = 0.166311; deltaR[26] = -0.0405203; deltaPhi[ 26] = -0.0145529; deltaZ[26] = 0;
// antTimeOffset[27] = 0.139405; deltaR[27] = -0.00401756; deltaPhi[ 27] = -0.00666063; deltaZ[27] = 0;
// antTimeOffset[28] = 0.104575; deltaR[28] = -0.0362955; deltaPhi[ 28] = -0.00372999; deltaZ[28] = 0;
// antTimeOffset[29] = 0.0345155; deltaR[29] = -0.00587152; deltaPhi[ 29] = 0.00197442; deltaZ[29] = 0;
// antTimeOffset[30] = 0.082859; deltaR[30] = -0.00611182; deltaPhi[ 30] = -0.000789595; deltaZ[30] = 0;
// antTimeOffset[31] = 0.078192; deltaR[31] = -0.00321244; deltaPhi[ 31] = 0.000188257; deltaZ[31] = 0;
// antTimeOffset[32] = 0.0143052; deltaR[32] = -0.0437687; deltaPhi[ 32] = -0.00289577; deltaZ[32] = 0;
// antTimeOffset[33] = 0.0671887; deltaR[33] = -0.0643475; deltaPhi[ 33] = -0.0203117; deltaZ[33] = 0;
// antTimeOffset[34] = 0.0819606; deltaR[34] = -0.0804245; deltaPhi[ 34] = -0.00503387; deltaZ[34] = 0;
// antTimeOffset[35] = -0.0659071; deltaR[35] = -0.0112675; deltaPhi[ 35] = -0.000220575; deltaZ[35] = 0;
// antTimeOffset[36] = -0.212624; deltaR[36] = 0.0337428; deltaPhi[ 36] = -0.00416114; deltaZ[36] = 0;
// antTimeOffset[37] = -0.00431668; deltaR[37] = -0.0525977; deltaPhi[ 37] = -0.0223176; deltaZ[37] = 0;
// antTimeOffset[38] = 0.0899397; deltaR[38] = -0.101587; deltaPhi[ 38] = 0.0058874; deltaZ[38] = 0;
// antTimeOffset[39] = 0.0294537; deltaR[39] = -0.0401037; deltaPhi[ 39] = 0.00899651; deltaZ[39] = 0;

// antTimeOffset[0] = 0.16044; deltaR[0] = -0.0618352; deltaPhi[ 0] = -0.00280991; deltaZ[0] = 0.0287513;
// antTimeOffset[1] = 0.290365; deltaR[1] = -0.110188; deltaPhi[ 1] = -0.00372876; deltaZ[1] = 0.0170052;
// antTimeOffset[2] = 0.349607; deltaR[2] = -0.140583; deltaPhi[ 2] = -0.024769; deltaZ[2] = 0.0623979;
// antTimeOffset[3] = 0.033917; deltaR[3] = -0.0260668; deltaPhi[ 3] = 0.00113633; deltaZ[3] = 0.0423352;
// antTimeOffset[4] = -0.074643; deltaR[4] = -0.0356161; deltaPhi[ 4] = -0.0141947; deltaZ[4] = -0.00394529;
// antTimeOffset[5] = 0.156259; deltaR[5] = -0.095648; deltaPhi[ 5] = 0.0117896; deltaZ[5] = 0.0105241;
// antTimeOffset[6] = 0.212836; deltaR[6] = -0.0907411; deltaPhi[ 6] = -0.0122972; deltaZ[6] = 0.0516748;
// antTimeOffset[7] = 0.541271; deltaR[7] = -0.157934; deltaPhi[ 7] = -0.00368077; deltaZ[7] = 0.130464;
// antTimeOffset[8] = -0.128575; deltaR[8] = 0.0582477; deltaPhi[ 8] = -0.0111813; deltaZ[8] = -0.035067;
// antTimeOffset[9] = -0.122485; deltaR[9] = -0.0119241; deltaPhi[ 9] = -0.00668961; deltaZ[9] = -0.0579302;
// antTimeOffset[10] = -0.263916; deltaR[10] = 0.0811694; deltaPhi[ 10] = -0.00717066; deltaZ[10] = -0.0897962;
// antTimeOffset[11] = -0.229033; deltaR[11] = 0.00897635; deltaPhi[ 11] = -5.65709e-05; deltaZ[11] = -0.0366053;
// antTimeOffset[12] = -0.222093; deltaR[12] = 0.0800239; deltaPhi[ 12] = -0.00901442; deltaZ[12] = -0.0750748;
// antTimeOffset[13] = -0.68459; deltaR[13] = 0.108285; deltaPhi[ 13] = -0.0148342; deltaZ[13] = -0.133195;
// antTimeOffset[14] = -0.516002; deltaR[14] = 0.107793; deltaPhi[ 14] = -0.00133016; deltaZ[14] = -0.143198;
// antTimeOffset[15] = -0.749822; deltaR[15] = 0.145538; deltaPhi[ 15] = -0.00448024; deltaZ[15] = -0.118927;
// antTimeOffset[16] = 0.0961751; deltaR[16] = -0.0564912; deltaPhi[ 16] = -0.00220892; deltaZ[16] = 0.0532255;
// antTimeOffset[17] = -0.166225; deltaR[17] = 0.00298098; deltaPhi[ 17] = 0.00128378; deltaZ[17] = -0.0263557;
// antTimeOffset[18] = 0.178941; deltaR[18] = -0.0972098; deltaPhi[ 18] = -0.00837075; deltaZ[18] = 0.0461125;
// antTimeOffset[19] = 0.171386; deltaR[19] = -0.0643506; deltaPhi[ 19] = -0.00831496; deltaZ[19] = 0.0592813;
// antTimeOffset[20] = -0.0204964; deltaR[20] = -0.0803763; deltaPhi[ 20] = -0.00428327; deltaZ[20] = -0.00602899;
// antTimeOffset[21] = 0.160323; deltaR[21] = -0.101956; deltaPhi[ 21] = -0.0151547; deltaZ[21] = 0.0514858;
// antTimeOffset[22] = 0.053953; deltaR[22] = -0.0773129; deltaPhi[ 22] = -0.00414444; deltaZ[22] = 0.0562807;
// antTimeOffset[23] = 0.00325878; deltaR[23] = -0.0378326; deltaPhi[ 23] = 0.000528826; deltaZ[23] = 0.0517299;
// antTimeOffset[24] = -0.0670197; deltaR[24] = -0.0225655; deltaPhi[ 24] = -0.00391848; deltaZ[24] = 0.0263552;
// antTimeOffset[25] = 0.0035097; deltaR[25] = -0.0427118; deltaPhi[ 25] = -0.00285503; deltaZ[25] = 0.0211308;
// antTimeOffset[26] = -0.0337735; deltaR[26] = -0.062921; deltaPhi[ 26] = -0.00747562; deltaZ[26] = -0.0138518;
// antTimeOffset[27] = -0.110361; deltaR[27] = -0.0144298; deltaPhi[ 27] = -0.0028118; deltaZ[27] = -0.0318333;
// antTimeOffset[28] = -0.227861; deltaR[28] = -0.0389804; deltaPhi[ 28] = -0.0012114; deltaZ[28] = -0.0674168;
// antTimeOffset[29] = -0.081439; deltaR[29] = -0.0449712; deltaPhi[ 29] = -0.00826954; deltaZ[29] = 0.00548985;
// antTimeOffset[30] = 0.17098; deltaR[30] = -0.0938958; deltaPhi[ 30] = -0.00969672; deltaZ[30] = 0.0597671;
// antTimeOffset[31] = 0.0258336; deltaR[31] = -0.052614; deltaPhi[ 31] = -0.000359792; deltaZ[31] = 0.044548;
// antTimeOffset[32] = -0.425254; deltaR[32] = -0.0372025; deltaPhi[ 32] = -0.00536103; deltaZ[32] = -0.0978299;
// antTimeOffset[33] = -0.440596; deltaR[33] = -0.0534945; deltaPhi[ 33] = -0.0144455; deltaZ[33] = -0.136331;
// antTimeOffset[34] = -0.256048; deltaR[34] = -0.0884831; deltaPhi[ 34] = -0.00580407; deltaZ[34] = -0.0676089;
// antTimeOffset[35] = -0.601107; deltaR[35] = 0.0224642; deltaPhi[ 35] = -0.00283517; deltaZ[35] = -0.115955;
// antTimeOffset[36] = -0.472965; deltaR[36] = 0.0023781; deltaPhi[ 36] = -0.00319073; deltaZ[36] = -0.049661;
// antTimeOffset[37] = -0.498374; deltaR[37] = -0.0241109; deltaPhi[ 37] = -0.0168551; deltaZ[37] = -0.102832;
// antTimeOffset[38] = -0.466373; deltaR[38] = -0.0623884; deltaPhi[ 38] = 0.00665517; deltaZ[38] = -0.131102;
// antTimeOffset[39] = -0.382128; deltaR[39] = -0.040161; deltaPhi[ 39] = 0.00508923; deltaZ[39] = -0.100047;

   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
   fGeomTool->useKurtAnitaIINumbers(1);
   char eventName[FILENAME_MAX];
   char headerName[FILENAME_MAX];
   char gpsName[FILENAME_MAX];
   char corrName[FILENAME_MAX];
   char outName[FILENAME_MAX];

   //   sprintf(baseDir,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/");
   sprintf(eventName,"%s/run%d/eventFile%d.root",baseDir,run,run);
   sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
   sprintf(gpsName,"%s/run%d/gpsEvent%d.root",baseDir,run,run);
   sprintf(corrName,"%s/corRun%d.root",corTreeDir,run);
   sprintf(outName,"%s/deltaTFile%d.root",outputDir,run);

   
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
   UInt_t expTaylorTime;
   TFile *fpCor = new TFile(corrName);
   TTree *corTree = (TTree*) fpCor->Get("corTree");
   corTree->SetBranchAddress("cor",&corSum);
   corTree->SetBranchAddress("labChip",&labChip);
   corTree->SetBranchAddress("expTaylorTime",&expTaylorTime);

   Long64_t numEntries=corTree->GetEntries();
   int counter=0;

   TFile *fpOut = new TFile(outName,"RECREATE");

   Long64_t entry=0;
   UInt_t eventNumber, triggerTime, triggerTimeNs;
   Int_t firstAnt,secondAnt,maxAnt,corInd;
   Double_t deltaT,deltaTExpected;
   Double_t phiWave, phiMaxAnt;
   Double_t thetaWave;
   Double_t corPeak, corRMS;
   Double_t balloonLat, balloonLon, balloonAlt;
   Double_t heading,pitch,roll;
   //   Double_t deltaZ;
   //   Double_t deltaR;
   Double_t meanPhiAntPair;
   Double_t deltaPhiAntPair;
   


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
   deltaTTree->Branch("phiMaxAnt",&phiMaxAnt,"phiMaxAnt/D");
   deltaTTree->Branch("phiWave",&phiWave,"phiWave/D");
   deltaTTree->Branch("thetaWave",&thetaWave,"thetaWave/D");
   deltaTTree->Branch("eventNumber",&eventNumber,"eventNumber/i");
   deltaTTree->Branch("triggerTime",&triggerTime,"triggerTime/i");
   deltaTTree->Branch("triggerTimeNs",&triggerTimeNs,"triggerTimeNs/i");
   deltaTTree->Branch("corInd",&corInd,"corInd/I");
   deltaTTree->Branch("balloonLat",&balloonLat,"balloonLat/D");
   deltaTTree->Branch("balloonLon",&balloonLon,"balloonLon/D");
   deltaTTree->Branch("balloonAlt",&balloonAlt,"balloonAlt/D");
   deltaTTree->Branch("heading",&heading,"heading/D");
   deltaTTree->Branch("pitch",&pitch,"pitch/D");
   deltaTTree->Branch("roll",&roll,"roll/D");
   //   deltaTTree->Branch("deltaZ",&deltaZ,"deltaZ/D");
   //   deltaTTree->Branch("deltaR",&deltaR,"deltaR/D");
   deltaTTree->Branch("meanPhiAntPair",&meanPhiAntPair,"meanPhiAntPair/D");
   deltaTTree->Branch("deltaPhiAntPair",&deltaPhiAntPair,"deltaPhiAntPair/D");
   deltaTTree->Branch("expTaylorTime",&expTaylorTime,"expTaylorTime/i");

  // Double_t thetaWave;

   for(entry=0;entry<numEntries;entry++) { 
  

      corTree->GetEntry(entry);
      Long64_t headEntry=headTree->GetEntryNumberWithIndex(corSum->eventNumber);
      if(headEntry<0) 
	continue;
      headTree->GetEntry(headEntry);
     
      // if(header->triggerTimeNs*1e-9< 0.097 || header->triggerTimeNs*1e-9>0.1)

      triggerTimeNs=header->triggerTimeNs;
      triggerTime=header->triggerTime;
      eventNumber=header->eventNumber;

      adu5PatTree->GetEntry(headEntry);
      
      
   // PrettyAnitaEvent realEvent(event,WaveCalType::kVTFullAGCrossCorClock,header);
      balloonLat=pat->latitude;
      balloonLon=pat->longitude;
      balloonAlt=pat->altitude;
      heading=pat->heading;
      pitch=pat->pitch;
      roll=pat->roll;

	//Simon 30/04/09 numbers     
      //      pat->pitch=0.64;
      //      pat->roll=0.14;


	//Simon 02/05/09 numbers     
      //      pat->pitch=0.76;
      //      pat->roll=0.13;

      //Test heading offset
      //      pat->heading+=0.19;
      //      if(pat->heading>=360) pat->heading-=360;
      //      if(pat->heading<0) pat->heading+=360;

     // 

      //      pat->pitch=0.5;
      //      pat->roll=-0.1;


      //Kurt numbers 
      pat->pitch=-0.29;
      pat->roll=0.89;
//       //Test heading offset
//      pat->heading-=0.32;
      if(pat->heading>=360) pat->heading-=360;
      if(pat->heading<0) pat->heading+=360;

      UsefulAdu5Pat usefulPat(pat);
     
      for(corInd=0;corInd<35;corInd++) {

	 firstAnt=corSum->firstAnt[corInd];
	 secondAnt=corSum->secondAnt[corInd];
	 deltaT=corSum->maxCorTimes[corInd];
	 maxAnt=corSum->centreAntenna;
	 phiMaxAnt=fGeomTool->getAntPhiPositionRelToAftFore(corSum->centreAntenna,AnitaPol::kHorizontal)*TMath::RadToDeg();
	 
	 //Default values
	 deltaTExpected=usefulPat.getDeltaTTaylor(corSum->firstAnt[corInd],corSum->secondAnt[corInd],AnitaPol::kHorizontal);
	 deltaTExpected=usefulPat.getDeltaTWillySeavey(corSum->firstAnt[corInd],corSum->secondAnt[corInd],AnitaPol::kHorizontal);

	 //These two lines are for use with Simon's array of numbers
	 //	 deltaTExpected=usefulPat.getDeltaTTaylorOpt(corSum->firstAnt[corInd],corSum->secondAnt[corInd],deltaR,deltaZ,deltaPhi);
	 //	 deltaT+=(antTimeOffset[firstAnt]-antTimeOffset[secondAnt]);



	 phiWave=usefulPat.getPhiWave();
	 thetaWave=usefulPat.getThetaWave();
	 corPeak=corSum->maxCorVals[corInd];
	 corRMS=corSum->rmsCorVals[corInd];
	 
	 meanPhiAntPair=fGeomTool->getMeanAntPairPhiRelToAftFore(firstAnt,secondAnt,AnitaPol::kHorizontal);
	 deltaPhiAntPair=fGeomTool->getPhiDiff(phiWave,meanPhiAntPair);//,AnitaPol::kHorizontal);
	 
	 //	 std::cout << phiWave << "\t" << meanPhiAntPair << "\t" << deltaPhiAntPair << "\n";

	 //Convert to degrees
	 deltaPhiAntPair*=TMath::RadToDeg();
	 phiWave*=TMath::RadToDeg();
	 thetaWave*=TMath::RadToDeg();	       

	 //Actually fill the tree
	 deltaTTree->Fill();
      }


      counter++; 
      if(counter%100==0)
	 cerr << "*";
   }
   deltaTTree->AutoSave();
   fpOut->Close();
   //   histSimpleDtDiff->Draw();
}

