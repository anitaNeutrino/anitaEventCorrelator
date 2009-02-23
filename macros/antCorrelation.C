

void antCorrelation() 
{

  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  gSystem->AddIncludePath("-L${ANITA_UTIL_INSTALL_DIR}/lib");
  
  gSystem->Load("libMathMore.so");
  gSystem->Load("/sw/lib/libfftw3.so");
  gSystem->Load("libRootFftwWrapper.so");
  gSystem->Load("/Applications/root/lib/libMinuit.so");
  gSystem->Load("libAnitaEvent.so");
   gSystem->Load("libAnitaCorrelator.so");
   TChain *deltaTTree = new TChain("deltaTTree");
   deltaTTree->Add("../outFiles/deltaTFile18.root");
   
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();

   double deltaR[16] = {0};
   double deltaZ[16] = {0};
   double deltaPhi[16] = {0};
   double deltaHeading[1] = {0};

   for(int i = 0; i < 16; i++){
     deltaR[i] =  0 ;
   }

   deltaHeading[0]=0;

  // deltaPhi[0] = 0.0494191 ;
//  deltaPhi[1] = 0.00782821;
//  deltaPhi[2] = -0.00444051;
//  deltaPhi[3] = -0.0340429;
//  deltaPhi[4] = -0.011112;
//  deltaPhi[5] = -0.0169422;
//  deltaPhi[6] = 0.0385532;
//  deltaPhi[7] =  0.00413145;
//  deltaPhi[8] = -0.0148806;
//  deltaPhi[9] = 0.0266541;
//  deltaPhi[10] = 0.0329093;
//  deltaPhi[11] = -0.0140619;
//  deltaPhi[12] = 0.0501811;
//  deltaPhi[13] = -0.000886006;
//  deltaPhi[14] =  0.0243596;
//  deltaPhi[15] = 0.0168409;


gSystem->CompileMacro("correlationTreeLoopOpt.C","k");
// correlationTreeLoopOpt(18,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/","../../outFiles/","../../outFiles/",deltaR,deltaZ,deltaPhi,deltaHeading);


correlationTreeLoopOpt(18,"/unix/anita3/flight0809/root/","../outFiles/","../outFiles/",deltaR,deltaZ,deltaPhi,deltaHeading);

   char plotCond[180];
   char plotTitle[180];
   char histName[180];
   char histNameOpt[180];
   TF1 *fitty = new TF1("fitty","gaus",-3,3);
   
   TH1F *histAttenDiff[32][4];
   TH1F *histAttenDiffOpt[32][4];
   int middleAnt; 
       int leftAnt,rightAnt;

   for(int chip=0; chip<4;chip++){
     for(int atten = 0; atten <32; atten++){
       middleAnt = atten; 
       fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 
       sprintf(plotTitle,"ANT %d and ANT %d",middleAnt, rightAnt);
	 sprintf(histName,"histAttenDiff%d_%d_%d",atten,atten+9,chip);
	 sprintf(histNameOpt,"histAttenDiffOpt%d_%d_%d",atten,atten+1,chip);
	 histAttenDiff[atten][chip]=new TH1F(histName,plotTitle,200,-1,1);
	 histAttenDiffOpt[atten][chip]=new TH1F(histNameOpt,plotTitle,200,-1,1);
     }
   }
   
   for(int chip=0; chip<1;chip++){
    for(int ant=0;ant<32;ant++) {
        middleAnt=ant; 
       //int leftAnt,rightAnt;
       fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 

//       //fitty->SetParLimits(2,0,1);
// 	 // sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",leftAnt,middleAnt,chip);	 
// 	 //sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",middleAnt,leftAnt,chip);	 
       //sprintf(plotCond,"((firstAnt==%d && secondAnt==%d) && labChip==%d) ",middleAnt,rightAnt,chip);	 
       sprintf(plotCond,"((firstAnt==%d && secondAnt==%d)) ",middleAnt,rightAnt);	 
 	 sprintf(plotTitle,"Ant %d -  Ant %d",middleAnt,rightAnt);
 	 sprintf(histName,"histDt_%d_%d",ant,chip);
 	 sprintf(histNameOpt,"histDtOpt_%d_%d",ant,chip);
	 if(chip==0){
	    std::cout << plotTitle << std::endl;
 	 }
		 TH1F *histDt11 = new TH1F(histName,plotTitle,200,-1,1);
		 TH1F *histDt11Opt = new TH1F(histNameOpt,plotTitle,200,-1,1);
 	 deltaTTree->Project(histName,"deltaT-deltaTExpected",plotCond);
	 deltaTTree->Project(histNameOpt,"deltaT-deltaTExpectedOpt",plotCond);
 	 //deltaTTree->Project(histName,"deltaT",plotCond);
	 histAttenDiff[ant][chip]->Add(histDt11);
	 histAttenDiffOpt[ant][chip]->Add(histDt11Opt);
    }
   }


     ofstream Output("diffSurfFastClock.txt");
    TCanvas *canSurf= new TCanvas("canSurf","canSurf");
    


    canSurf->Divide(4,4);

    //TH1F *addedHist = new TH1F(histNameOpt,plotTitle,200,-1,1);
    Double_t histAdded = 0;
    Double_t histAddedRMS = 0;
    Double_t histAddedOpt = 0;
    Double_t histAddedRMSOpt = 0;

    for(int atten=0;atten<16;atten++) {
      for(int chip=0; chip<1;chip++){
	canSurf->cd(atten+1);
	histAttenDiff[atten][chip]->SetLineColor(chip+1);
	histAttenDiffOpt[atten][chip]->SetLineColor(chip+2);
	if(chip==0){
	  histAttenDiff[atten][chip]->Draw();
	   histAttenDiffOpt[atten][chip]->Draw("same");

	  histAdded = histAdded + histAttenDiff[atten][chip]->GetMean();
	  histAddedRMS = histAddedRMS + histAttenDiff[atten][chip]->GetRMS();
	  histAddedOpt = histAddedOpt + histAttenDiffOpt[atten][chip]->GetMean();
	  histAddedRMSOpt = histAddedRMSOpt + histAttenDiffOpt[atten][chip]->GetRMS();
	}
	else{
	  histAttenDiff[atten][chip]->Draw("same");
	  histAdded = histAdded + histAttenDiff[atten][chip]->GetMean();
	  histAddedRMS = histAddedRMS + histAttenDiff[atten][chip]->GetRMS();
	}      
      }
    }

    cout << histAdded << "  "<< histAddedRMS << endl;
    cout << histAddedOpt << "  "<< histAddedRMSOpt << endl;


//   TCanvas *canSurf= new TCanvas("canSurf","canSurf");
//     canSurf->Divide(1,2);



//  	canSurf->cd(1);
//  	histAttenDiff[0][0]->SetLineColor(1);
//  	histAttenDiffOpt[0][0]->SetLineColor(2);
 	
	//histAttenDiff[0][0]->Draw();
	//canSurf->cd(2); 	
	//histAttenDiffOpt[0][0]->Draw();
    
}
