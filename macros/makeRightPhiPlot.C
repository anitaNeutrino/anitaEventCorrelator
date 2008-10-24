

void makeRightPhiPlot() 
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
   deltaTTree->Add("deltaTFileFast*.root");
   
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();



   char plotCond[180];
   char plotTitle[180];
   char histName[180];
   TF1 *fitty = new TF1("fitty","gaus",-3,3);

   TH1F *histMeanDiff = new TH1F("histMeanDiff","histMeanDiff",100,-3,3);

   for(int ant=0;ant<16;ant++) {
      TCanvas *can = new TCanvas();//"can","can");
      can->Divide(2,2);
      int middleAnt=ant; 
      int leftAnt,rightAnt;
      fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 
      cout << leftAnt << "\t" << middleAnt << "\t" << rightAnt << "\n";
      for(int chip=0;chip<4;chip++) {
	 fitty->SetParLimits(2,0,1);

	 can->cd(chip+1);
	 //	 sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",rightAnt,middleAnt,chip);	 
	 sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",middleAnt,rightAnt,chip);	 
	 //	 sprintf(plotTitle,"Ant %d -  Ant %d (Chip %d)",rightAnt,middleAnt,chip);
	 sprintf(plotTitle,"Ant %d -  Ant %d (Chip %d)",middleAnt,rightAnt,chip);
	 sprintf(histName,"histDt_%d_%d",ant,chip);
	 TH1F *histDt11 = new TH1F(histName,plotTitle,100,-3,3);

	 deltaTTree->Project(histName,"deltaT-deltaTExpected",plotCond);
	 //	 cout << plotCond << endl;
	 histDt11->Draw();
	 fitty->SetRange(histDt11->GetMean()-0.5,histDt11->GetMean()+0.5);
	 fitty->SetParameters(10,histDt11->GetMean(),0.1);
	 histDt11->Fit("fitty","QR");
	 Int_t numUnder= histDt11->GetBinContent(0);
	 Int_t numOver =histDt11->GetBinContent(1+histDt11->GetNbinsX());
	 cout << rightAnt << "\t" << middleAnt << "\t" << chip << "\t" << histDt11->GetEntries() << "\t"
	      << (histDt11->GetEntries()-(numOver+numUnder)) << "\t" << histDt11->GetMean() << "\t" << histDt11->GetRMS() << "\t"
	      << fitty->GetParameter(1) << "\t" << fitty->GetParError(1) << "\t" << fitty->GetParameter(2) << "\t"
	      << fitty->GetParError(2)  << "\t" << fitty->GetChisquare() << "\t" << fitty->GetNDF() << "\n";
	 histMeanDiff->Fill(fitty->GetParameter(1));
	 
      }
	 
   }
   //   TCanvas *canSum = new TCanvas();
   //   histMeanDiff->Draw();
      
}
