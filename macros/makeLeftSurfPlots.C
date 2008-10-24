

void makeLeftSurfPlots() 
{
   gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
   gSystem->AddIncludePath("-I${PLOTTER_DIR}");
   //  cout << gSystem->GetIncludePath() <<endl;
   
   gSystem->Load("libMathMore.so");
   gSystem->Load("/usr/lib64/libfftw3.so");
   gSystem->Load("libAnitaEvent.so");
   gSystem->Load("libAnitaCorrelator.so");
   TChain *deltaTTree = new TChain("deltaTTree");
   //   TFile *fp = new TFile("deltaTFile1027.root");
   //   TTree *deltaTTree = (TTree*) fp->Get("deltaTTree");
   deltaTTree->Add("deltaTFileFastClock*.root");
   
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();



   char plotCond[180];
   char plotTitle[180];
   char histName[180];
   TF1 *fitty = new TF1("fitty","gaus",-3,3);
   
   TH1F *histSurfDiff[8][4];
   for(int surf=0;surf<8;surf++) {
      for(int chip=0;chip<4;chip++) {
	 sprintf(plotTitle,"SURF %d - SURF %d (Chip %d)",surf,(surf+1)%8,chip);
	 sprintf(histName,"histSurfDiff%d_%d",surf,chip);
	 histSurfDiff[surf][chip]=new TH1F(histName,plotTitle,100,-3,3);
      }
   }
   
   for(int ant=0;ant<32;ant++) {
      int middleAnt=ant; 
      int leftAnt,rightAnt;
      fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 
      Int_t ci=fGeomTool->getChanIndexFromAntPol(leftAnt,AnitaPol::kVertical);
      Int_t firstSurf=ci/9;

      for(int chip=0;chip<4;chip++) {
	 fitty->SetParLimits(2,0,1);
	 sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",leftAnt,middleAnt,chip);	 
	 //	 sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",middleAnt,leftAnt,chip);	 
	 sprintf(plotTitle,"Ant %d -  Ant %d (Chip %d)",leftAnt,middleAnt,chip);
	 sprintf(histName,"histDt_%d_%d",ant,chip);
	 TH1F *histDt11 = new TH1F(histName,plotTitle,100,-3,3);
	 deltaTTree->Project(histName,"deltaT-deltaTExpected",plotCond);
	 histSurfDiff[firstSurf][chip]->Add(histDt11);
      }	 
   }
   //   ofstream Output("diffSurfFastClock.txt");
   TCanvas *canSurf= new TCanvas("canSurf","canSurf");
   canSurf->Divide(2,4);
   for(int surf=0;surf<8;surf++) {
      canSurf->cd(surf+1);
      for(int chip=0;chip<4;chip++) {
	 histSurfDiff[surf][chip]->SetLineColor(getNiceColour(chip));
	 if(chip==0)
	    histSurfDiff[surf][chip]->Draw();
	 else
	    histSurfDiff[surf][chip]->Draw("same");
	    
// 	 histSurfDiff[surf][chip]->Fit("fitty","QR");
// 	 Output << surf << "\t" << (surf+1)%8 << "\t" << chip << "\t" << fitty->GetParameter(1) << "\t" << fitty->GetParError(1) << "\n";
      }
   }
   
}
