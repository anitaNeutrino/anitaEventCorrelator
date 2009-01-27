

void attenCorrelation() 
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
   deltaTTree->Add("../../outFiles/deltaTFileSlowClock18.root");
   
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();

   char plotCond[180];
   char plotTitle[180];
   char histName[180];
   TF1 *fitty = new TF1("fitty","gaus",-3,3);
   
   TH1F *histAttenDiff[32];
   
   for(int atten = 0; atten <32; atten++){
	 sprintf(plotTitle,"ATTEN %d - ATTEN %d",atten,atten+1);
	 sprintf(histName,"histAttenDiff%d_%d",atten,atten+1);
	 histAttenDiff[atten]=new TH1F(histName,plotTitle,500,-3,3);
      }
   
   
    for(int ant=0;ant<32;ant++) {
       int middleAnt=ant; 
       int leftAnt,rightAnt;
       fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 

//       //fitty->SetParLimits(2,0,1);
// 	 // sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",leftAnt,middleAnt,chip);	 
// 	 //sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",middleAnt,leftAnt,chip);	 
 	 sprintf(plotCond,"((firstAnt==%d && secondAnt==%d) && labChip==3) ",middleAnt,rightAnt);	 
 	 sprintf(plotTitle,"Ant %d -  Ant %d",leftAnt,middleAnt);
 	 sprintf(histName,"histDt_%d",ant);
	 //	 std::cout << plotTitle << std::endl;
 	 TH1F *histDt11 = new TH1F(histName,plotTitle,500,-3,3);
 	 deltaTTree->Project(histName,"deltaT-deltaTExpected",plotCond);
 	 //deltaTTree->Project(histName,"deltaT",plotCond);
	 histAttenDiff[ant]->Add(histDt11);
    }



   //   ofstream Output("diffSurfFastClock.txt");
    TCanvas *canSurf= new TCanvas("canSurf","canSurf");
    canSurf->Divide(4,8);
    for(int atten=0;atten<32;atten++) {
      canSurf->cd(atten+1);
      histAttenDiff[atten]->SetLineColor(1);
      histAttenDiff[atten]->Draw();
	    
    }
   
   
}
