

Double_t getDelay(Double_t antPhi) {
  return 1e-8*TMath::Abs(antPhi*antPhi*antPhi*antPhi); //
}

Double_t getLeftPhi(Double_t midPhi)
{return (midPhi+11.25);}
Double_t getRightPhi(Double_t midPhi)
{return (midPhi-11.25);}

void plotExampleSShape() {
  plotAllLowerAnts();
  //  justPlotFunction();

}

Double_t singleDelayWithGeom(Double_t *x, Double_t *par) {
  Double_t antPhi=x[0];
  Double_t antDelay=singleDelayFunc(x,par);
  Double_t geomDelay=0.2e9*(1-TMath::Cos(TMath::DegToRad()*antPhi))/299792458.;
  return antDelay+geomDelay;
}

Double_t singleDelayFunc(Double_t *x, Double_t *par) {
  Double_t antPhi=x[0];
  Double_t antDelay=par[1]*antPhi + par[2]*antPhi*antPhi + par[3]*antPhi*antPhi*antPhi + par[4]*antPhi*antPhi*antPhi*antPhi;
  return antDelay;
}

Double_t groupDelayFunc(Double_t *x, Double_t *par) {
  Double_t leftPhi=x[0]+11.25;
  Double_t rightPhi=x[0]-11.25;
  
  //leftPhi*=-1;
  //  rightPhi*=-1;
  Double_t leftDelay=singleDelayFunc(&leftPhi,par);
  Double_t rightDelay=singleDelayFunc(&rightPhi,par);

  return (leftDelay-rightDelay);
}

void plotAllLowerAnts() 
{
  gStyle->SetFitFormat("6.3e");
  gSystem->Load("libAnitaEvent.so");
  TChain *deltaChain = new TChain("deltaTTree");
  deltaChain->Add("/home/rjn/anita/data/deltaTTrees/justTaylorSvn220509/deltaTFile*.root");
  TF1 *fitDelay = new TF1("fitDelay",groupDelayFunc,-50,50,5);
  char plotCond[180];  

  TCanvas *canData2 = new TCanvas("canData2","canData2",600,400); 
  sprintf(plotCond,"(firstAnt>=16 && firstAnt<32 && secondAnt==(firstAnt+1) && corPeak/corRMS>8) && corRMS>40");
  deltaChain->Draw("deltaT-deltaTExpected:deltaPhiAntPair",plotCond,"prof");
  {
    TH1F *htemp = (TH1F*)gPad->FindObject("htemp");
    if(htemp) {
      htemp->GetXaxis()->SetRangeUser(-50,50);
      htemp->GetYaxis()->SetRangeUser(-0.1,0.1);
      fitDelay->SetParameters(0,0,0,0,1e-8);
      fitDelay->SetParLimits(0,-1,-1);
      fitDelay->SetParLimits(1,-1,-1);
      //      fitDelay->SetParLimits(2,-1,-1);
      fitDelay->SetParLimits(3,-1,-1);
      htemp->Fit("fitDelay","R");
      htemp->SetTitle("Phi S-Shape");
      htemp->GetXaxis()->SetTitle("Phi Wave - Mid Ant (Degrees)");
      htemp->GetYaxis()->SetTitle("#Delta#Deltat (ns)");
    }
  } 
  
  TCanvas *canDelay = new TCanvas("canDelay","canDelay",600,600);  
  canDelay->Divide(1,2);
  {
    canDelay->cd(1);
    TF1 *antDelay = new TF1("antDelay",singleDelayFunc,-50,50,5);
    antDelay->SetParameters(fitDelay->GetParameters());
    antDelay->Draw();    
    antDelay->GetXaxis()->SetTitle("Phi from Boresight");
    antDelay->GetYaxis()->SetTitle("Delay (ns)");
  }
  {
    canDelay->cd(2);
    TF1 *antDelay = new TF1("antDelayGeom",singleDelayWithGeom,-50,50,5);
    antDelay->SetParameters(fitDelay->GetParameters());
    antDelay->Draw();    
    antDelay->GetXaxis()->SetTitle("Phi from Boresight");
    antDelay->GetYaxis()->SetTitle("Delay + Geom (ns)");
  }
  

//   TCanvas *canData = new TCanvas("canData","canData",800,800);  
//   canData->Divide(4,4);
 
//   for(int ant=16;ant<32;ant++) {
//     canData->cd(ant-15);
//     fitDelay->SetParameters(0,0,0,0,1e-8);
//     fitDelay->SetParLimits(0,-1,-1);
//     //    sprintf(plotCond,"firstAnt>=16 && firstAnt<32 && secondAnt==(firstAnt+1) && corPeak/corRMS>8 ");
//     sprintf(plotCond,"firstAnt==%d && firstAnt<32 && secondAnt==(firstAnt+1) && corPeak/corRMS>8 ",ant);
//     if(ant==31)
//       sprintf(plotCond,"firstAnt==31 && secondAnt==16 && corPeak/corRMS>8");
	    
//     deltaChain->Draw("deltaT-deltaTExpected:deltaPhiAntPair",plotCond,"prof");
//     {
//       TH1F *htemp = (TH1F*)gPad->FindObject("htemp");
//       if(htemp) {
// 	htemp->GetXaxis()->SetRangeUser(-50,50);
// 	htemp->GetYaxis()->SetRangeUser(-0.1,0.1);
// 	htemp->Fit("fitDelay","R");
// 	htemp->SetTitle("Phi S-Shape");
// 	htemp->GetXaxis()->SetTitle("Phi Wave - Mid Ant (Degrees)");
// 	htemp->GetYaxis()->SetTitle("#Delta#Deltat (ns)");
//       }
//     }
//     //  fitDelay->Draw("same");
//   }
}


void justPlotFunction() {
  //Here we are looking at left -right 
  Double_t midPhi[1000];
  Double_t leftPhi[1000];
  Double_t rightPhi[1000];
  Double_t leftTime[1000];
  Double_t rightTime[1000];
  Double_t diffTime[1000];
    
  Int_t numPoints=0;
  for(Double_t phi=-50;phi<50;phi+=0.1) {
    midPhi[numPoints]=phi;
    leftPhi[numPoints]=getLeftPhi(phi);
    rightPhi[numPoints]=getRightPhi(phi);
    leftTime[numPoints]=getDelay(leftPhi[numPoints]);
    rightTime[numPoints]=getDelay(rightPhi[numPoints]);
    diffTime[numPoints]=leftTime[numPoints]-rightTime[numPoints];
    numPoints++;
  }
  TCanvas *canPhi = new TCanvas("canPhi","canPhi");
  TMultiGraph *mgPhi = new TMultiGraph();
  TGraph *grPhiLeft = new TGraph(numPoints,midPhi,leftPhi);
  grPhiLeft->SetMarkerColor(getNiceColour(0));
  grPhiLeft->SetMarkerStyle(getMarker(0));
  mgPhi->Add(grPhiLeft,"p");
  TGraph *grPhiRight = new TGraph(numPoints,midPhi,rightPhi);
  grPhiRight->SetMarkerColor(getNiceColour(1));
  grPhiRight->SetMarkerStyle(getMarker(1));
  mgPhi->Add(grPhiRight,"p");
  mgPhi->Draw("ap");

  TCanvas *canTime = new TCanvas("canTime","canTime",800,800);
  canTime->Divide(1,2);
  canTime->cd(1);
  TMultiGraph *mgTime = new TMultiGraph();
  TGraph *grTimeLeft = new TGraph(numPoints,midPhi,leftTime);
  grTimeLeft->SetMarkerColor(getNiceColour(0));
  grTimeLeft->SetMarkerStyle(getMarker(0));
  mgTime->Add(grTimeLeft,"p");
  TGraph *grTimeRight = new TGraph(numPoints,midPhi,rightTime);
  grTimeRight->SetMarkerColor(getNiceColour(1));
  grTimeRight->SetMarkerStyle(getMarker(1));
  mgTime->Add(grTimeRight,"p");
  mgTime->Draw("ap");
  canTime->cd(2);
  TGraph *grDiff = new TGraph(numPoints,midPhi,diffTime);
  grDiff->SetMarkerColor(getNiceColour(2));
  grDiff->SetMarkerStyle(getMarker(2));
  grDiff->Draw("ap");
}
