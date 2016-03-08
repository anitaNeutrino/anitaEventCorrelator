#define groupDelay_cxx
#include "groupDelay.h"
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>


//Here are the optimal Kurt Time Offsets
Double_t antTimeDiff[40]={0};
//Double_t antTimeDiff[40]={0.0933231,0.0983801,0.0319324,0.0943634,0.0244397,0.0955691,0.0414201,0.0661128,0.231548,0.0642528,0.267108,0.0407496,0.303581,0.0138057,0.166801,-0.00138648,0.0124803,0.0152421,-0.0359487,0.0281984,-0.091734,-0.0223349,-0.0412189,0.0384707,0.0325476,0.0564723,0.0198542,0.104875,-0.0662191,-0.0474661,-0.00594653,0.00272739,-0.269686,-0.245749,-0.238552,-0.209472,-0.235972,-0.233097,-0.272636,-0.214836};


Double_t getDelay(Double_t antPhi) {
  return 1e-8*TMath::Abs(antPhi*antPhi*antPhi*antPhi); //
}

Double_t getLeftPhi(Double_t midPhi)
{return (midPhi+11.25);}
Double_t getRightPhi(Double_t midPhi)
{return (midPhi-11.25);}

Double_t singleDelayFunc(Double_t *x, Double_t *par) {
  Double_t antPhi=x[0];
  Double_t antDelay=par[1]*antPhi + par[2]*antPhi*antPhi + par[3]*antPhi*antPhi*antPhi + par[4]*antPhi*antPhi*antPhi*antPhi;
  return antDelay;
}

Double_t singleDelayWithGeom(Double_t *x, Double_t *par) {
  Double_t antPhi=x[0];
  Double_t antDelay=singleDelayFunc(x,par);
  Double_t geomDelay=0.3e9*(1-TMath::Cos(TMath::DegToRad()*antPhi))/299792458.;
  return antDelay+geomDelay;
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


void groupDelay::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TH2F *histDtPhi = new TH2F("histDtPhi","#Delta#Deltat vs Mid Ant Phi (deg)",
			      100,-50,50,200,-1,1);

   char histName[180];
   char histTitle[180];
   TH2F *histDtPhiAnt[40];
   for(int ant=0;ant<40;ant++) {
     sprintf(histName,"histDtPhiAnt%d",ant);
     sprintf(histTitle,"#Delta#Deltat vs Mid Phi (deg)");
     histDtPhiAnt[ant]= new TH2F(histName,histTitle,100,-50,50,200,-1,1);
   }


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) < 0) continue;
      
      histDtPhi->Fill(deltaPhiAntPair,(deltaT+antTimeDiff[firstAnt]-antTimeDiff[secondAnt])-deltaTExpected);      
      histDtPhiAnt[firstAnt]->Fill(deltaPhiAntPair,(deltaT+antTimeDiff[firstAnt]-antTimeDiff[secondAnt])-deltaTExpected);      
   }


   TF1 *fitDelay = new TF1("fitDelay",groupDelayFunc,-50,50,5);
   fitDelay->SetParLimits(0,-1,-1);
   fitDelay->SetParLimits(1,-1,-1);
   fitDelay->SetParLimits(3,-1,-1);
   TCanvas *canData2 = new TCanvas("canData2","canData2",600,600); 
   canData2->Divide(1,2);
   canData2->cd(1);
   histDtPhi->Draw("colz");
   histDtPhi->GetXaxis()->SetTitle("Middle Phi (Degrees)");
   histDtPhi->GetYaxis()->SetTitle("#Delta#Deltat");
   histDtPhi->GetYaxis()->SetRangeUser(-0.2,0.2);
   canData2->cd(2);
   TProfile *histDtProf = histDtPhi->ProfileX();
   histDtProf->Draw();
   histDtProf->Fit("fitDelay");
   histDtProf->GetXaxis()->SetTitle("Middle Phi (Degrees)");
   histDtProf->GetYaxis()->SetTitle("#Delta#Deltat");
   histDtProf->GetYaxis()->SetRangeUser(-0.1,0.1);

   
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
    gStyle->SetFitFormat("6.3e");
    canDelay->cd(2);
    TF1 *antDelay = new TF1("antDelayGeom",singleDelayWithGeom,-50,50,5);
    antDelay->SetParameters(fitDelay->GetParameters());
    antDelay->Draw();    
    antDelay->GetXaxis()->SetTitle("Phi from Boresight");
    antDelay->GetYaxis()->SetTitle("Delay + Geom (ns)");
  }
  

  TCanvas *canAnt = new TCanvas("canAnt","canAnt",800,800);
  canAnt->Divide(4,10,0,0);
  for(int ant=0;ant<40;ant++) {
    canAnt->cd(ant+1);
    histDtPhiAnt[ant]->Draw("colz");
    histDtPhiAnt[ant]->GetXaxis()->SetTitle("Middle Phi (Degrees)");
    histDtPhiAnt[ant]->GetYaxis()->SetTitle("#Delta#Deltat");
    histDtPhiAnt[ant]->GetYaxis()->SetRangeUser(-0.2,0.2);
  }
 

  TCanvas *canProf = new TCanvas("canProf","canProf",800,800);
  canProf->Divide(4,10,0,0);
  for(int ant=0;ant<40;ant++) {
    canProf->cd(ant+1);
    TProfile *histDtProf = histDtPhiAnt[ant]->ProfileX();
    histDtProf->Draw();
    histDtProf->Fit("fitDelay");
    histDtProf->GetXaxis()->SetTitle("Middle Phi (Degrees)");
    histDtProf->GetYaxis()->SetTitle("#Delta#Deltat");
    histDtProf->GetYaxis()->SetRangeUser(-0.1,0.1);
  }
  
}
