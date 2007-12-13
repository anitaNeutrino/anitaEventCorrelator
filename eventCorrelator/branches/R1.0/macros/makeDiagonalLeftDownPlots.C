

void makeDiagonalLeftDownPlots() 
{
   gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
   gSystem->AddIncludePath("-I${PLOTTER_DIR}");
   //  cout << gSystem->GetIncludePath() <<endl;
   
   gSystem->Load("libMathMore.so");
   gSystem->Load("/usr/lib64/libfftw3.so");
   gSystem->Load("libAnitaEvent.so");
   gSystem->Load("libAnitaPlotter.so");
   TChain *deltaTTree = new TChain("deltaTTree");
   //   TFile *fp = new TFile("deltaTFile1027.root");
   //   TTree *deltaTTree = (TTree*) fp->Get("deltaTTree");
   deltaTTree->Add("deltaTFileClock*.root");
   
   AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();

   gStyle->SetStatX(1);
   gStyle->SetStatY(1);
   gStyle->SetStatW(0.3);
   gStyle->SetStatH(0.3);
   

   char plotCond[180];
   char plotTitle[180];
   char histName[180];
   TF1 *fitty = new TF1("fitty","gaus",-3,3);

   TH1F *histMeanDiff = new TH1F("histMeanDiff","histMeanDiff",1000,-0.5,0.5);
   TH1F *histChiSq = new TH1F("histChiSq","histChiSq",100,0,3);
   TH1F *histSigma = new TH1F("histSigma","histSigma",100,0,10);
   TH1F *histConstant = new TH1F("histConstant","histConstant",100,0,100);

   TCanvas *can = new TCanvas();//"can","can");
   can->Divide(8,8,0,0);

   ofstream Output("diffsDiagonalLeftDown.txt");


   for(int ant=0;ant<16;ant++) {     
      int middleAnt=ant; 
      int leftAnt,rightAnt;
      fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 
      //      cout << leftAnt << "\t" << middleAnt << "\t" << rightAnt << "\n";
      int leftDownAnt=fGeomTool->getAzimuthPartner(leftAnt);
      //      cout << leftDownAnt << "\n";


      for(int chip=0;chip<4;chip++) {
	 fitty->SetParLimits(2,0,1);	 
	 can->cd((4*ant)+chip+1);
	 //	 sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",leftDownAnt,middleAnt,chip);	 
	 sprintf(plotCond,"((firstAnt==%d && secondAnt==%d))  && labChip==%d && (corPeak/corRMS)>6",middleAnt,leftDownAnt,chip);	 
	 //	 sprintf(plotTitle,"Ant %d -  Ant %d (Chip %d)",leftDownAnt,middleAnt,chip);
	 sprintf(plotTitle,"Ant %d -  Ant %d (Chip %d)",middleAnt,leftDownAnt,chip);
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
// 	 cout << leftAnt << "\t" << middleAnt << "\t" << chip << "\t" << histDt11->GetEntries() << "\t"
// 	      << (histDt11->GetEntries()-(numOver+numUnder)) << "\t" << histDt11->GetMean() << "\t" << histDt11->GetRMS() << "\t"
// 	      << fitty->GetParameter(1) << "\t" << fitty->GetParError(1) << "\t" << fitty->GetParameter(2) << "\t"
// 	      << fitty->GetParError(2)  << "\t" << fitty->GetChisquare() << "\t" << fitty->GetNDF() << "\n";

	 Double_t constant=histDt11->GetMean();
	 Double_t error=0;
	 if(histDt11->GetEntries())
	    error=histDt11->GetRMS()/TMath::Sqrt(histDt11->GetEntries());
	 Int_t entries=(histDt11->GetEntries()-(numOver+numUnder));

	 if(TMath::Abs(fitty->GetParameter(1)-histDt11->GetMean())<0.01) {
	    constant=fitty->GetParameter(1);
	    error=fitty->GetParError(1);
	 }
	 Output << middleAnt << "\t" << leftDownAnt << "\t" << chip << "\t" << constant << "\t" << error << "\t"
		<< entries << "\n";

	 histMeanDiff->Fill(fitty->GetParameter(1)-histDt11->GetMean());
	 histChiSq->Fill(fitty->GetChisquare()/Double_t(fitty->GetNDF()));
	 histSigma->Fill(fitty->GetParameter(2));
	 histConstant->Fill(fitty->GetParameter(0));
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
