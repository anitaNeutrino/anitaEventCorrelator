

void makeSamePhiPlots() 
{
   gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
   //  cout << gSystem->GetIncludePath() <<endl;
   
   gSystem->Load("libMathMore.so");
   gSystem->Load("libMinuit.so");
   gSystem->Load("/sw/lib/libfftw3.so");
   gSystem->Load("libAnitaEvent.so");
   gSystem->Load("libAnitaPlotter.so");

   makeSamePhiPlotsFast("feedMinus0p2");
}


void makeSamePhiPlotsFast(char *dirName)
{
  char fileMatch[180];
  sprintf(fileMatch,"%s/deltaTFileFast*",dirName);
   TChain *deltaTTree = new TChain("deltaTTree");
   //   TFile *fp = new TFile("deltaTFile1027.root");
   //   TTree *deltaTTree = (TTree*) fp->Get("deltaTTree");
   deltaTTree->Add(fileMatch);
   
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();



   char plotCond[180];
   char plotTitle[180];
   char histName[180];
   TF1 *fitty = new TF1("fitty","gaus",-1,1);

   TH1F *histMeanDiff = new TH1F("histMeanDiff","histMeanDiff",100,-1,1);
   TH1F *histChiSq = new TH1F("histChiSq","histChiSq",100,0,3);
   TH1F *histSigma = new TH1F("histSigma","histSigma",100,0,10);
   TH1F *histConstant = new TH1F("histConstant","histConstant",100,0,100);

   ofstream Output("diffsDown.txt");

   for(int ant=0;ant<16;ant++) {
      TCanvas *can = new TCanvas();//"can","can");
      can->Divide(2,2);
      int topAnt=ant; 
      int bottomAnt=fGeomTool->getAzimuthPartner(topAnt); 
      for(int chip=0;chip<4;chip++) {
	 fitty->SetParameters(10,0,0.05);
	 //	 fitty->SetParLimits(2,0,0.1);

	 can->cd(chip+1);
	 sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",topAnt,bottomAnt,chip);	 
	 sprintf(plotTitle,"Ant %d -  Ant %d (Chip %d)",topAnt,bottomAnt,chip);
	 sprintf(histName,"histDt_%d_%d",ant,chip);
	 TH1F *histDt11 = new TH1F(histName,plotTitle,40,-0.5,0.5);

	 deltaTTree->Project(histName,"deltaT-deltaTExpected",plotCond);
	 //	 cout << plotCond << endl;
	 histDt11->Draw();
	 histDt11->Fit("fitty","Q");
	 Int_t numUnder= histDt11->GetBinContent(0);
	 Int_t numOver =histDt11->GetBinContent(1+histDt11->GetNbinsX());
// 	 cout << topAnt << "\t" << bottomAnt << "\t" << chip << "\t" << histDt11->GetEntries() << "\t"
// 	      << (histDt11->GetEntries()-(numOver+numUnder)) << "\t" << histDt11->GetMean() << "\t" << histDt11->GetRMS() << "\t"
// 	      << fitty->GetParameter(1) << "\t" << fitty->GetParError(1) << "\t" << fitty->GetParameter(2) << "\t"
// 	      << fitty->GetParError(2)  << "\t" << fitty->GetChisquare() << "\t" << fitty->GetNDF() << "\n";
	 histMeanDiff->Fill(fitty->GetParameter(1));
	 histChiSq->Fill(fitty->GetChisquare()/Double_t(fitty->GetNDF()));
	 histSigma->Fill(fitty->GetParameter(2));
	 histConstant->Fill(fitty->GetParameter(0));


	 Double_t constant=histDt11->GetMean();
	 Double_t error=0;
	 if(histDt11->GetEntries())
	    error=histDt11->GetRMS()/TMath::Sqrt(histDt11->GetEntries());
	 Int_t entries=(histDt11->GetEntries()-(numOver+numUnder));

	 if(TMath::Abs(fitty->GetParameter(1)-histDt11->GetMean())<0.01 || (entries>=20 && fitty->GetParError(1)<0.08)) {
	    constant=fitty->GetParameter(1);
	    error=fitty->GetParError(1);
	 }
	 else {
	    cout << topAnt << "\t" << bottomAnt << "\t" << chip << "\n";
	 }
	 Output << topAnt << "\t" << bottomAnt << "\t" << chip << "\t" << constant << "\t" << error << "\t"
		<< entries << "\n";
      }
	 
   }
   TCanvas *canSum = new TCanvas();
   canSum->Divide(2,2);
   canSum->cd(1);
   histChiSq->Draw();
   canSum->cd(2);
   histConstant->Draw();
   canSum->cd(3);
   histMeanDiff->Draw();
   canSum->cd(4);
   histSigma->Draw();    

      
}
