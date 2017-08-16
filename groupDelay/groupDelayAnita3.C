#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>


const int MAX_ANTENNAS = 48;

//Here are the optimal Kurt Time Offsets
Double_t antTimeDiff[MAX_ANTENNAS]={0};
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

Double_t singleDelayFuncMod(Double_t *x, Double_t *par) {
  Double_t antPhi=x[0];
  //  Double_t antDelay=par[1]*antPhi + par[2]*antPhi*antPhi + par[3]*antPhi*antPhi*antPhi + par[4]*antPhi*antPhi*antPhi*antPhi + par[5]*antPhi*antPhi*antPhi*antPhi*antPhi + par[6]*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi + par[7]*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi + par[8]*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi + par[9]*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi + par[10]*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi*antPhi;
  Double_t antDelay=par[1]*TMath::Power(antPhi, 1) +
    + par[2]*TMath::Power(antPhi, 2) 
    + par[3]*TMath::Power(antPhi, 3) 
    + par[4]*TMath::Power(antPhi, 4) 
    + par[5]*TMath::Power(antPhi, 5) 
    + par[6]*TMath::Power(antPhi, 6) 
    + par[7]*TMath::Power(antPhi, 7) 
    + par[8]*TMath::Power(antPhi, 8) 
    + par[9]*TMath::Power(antPhi, 9) 
    + par[10]*TMath::Power(antPhi, 10) 
    + par[11]*TMath::Power(antPhi, 11) 
    + par[12]*TMath::Power(antPhi, 12) 
    // + par[13]*TMath::Power(antPhi, 13) 
    // + par[14]*TMath::Power(antPhi, 14) 
    // + par[15]*TMath::Power(antPhi, 15) 
    // + par[16]*TMath::Power(antPhi, 16) 
    // + par[17]*TMath::Power(antPhi, 17) 
    // + par[18]*TMath::Power(antPhi, 18) 
    // + par[19]*TMath::Power(antPhi, 19) 
    // + par[20]*TMath::Power(antPhi, 20) 
    ;

  return antDelay;
}


Double_t singleDelayWithGeom(Double_t *x, Double_t *par) {
  Double_t antPhi=x[0];
  Double_t antDelay=singleDelayFunc(x,par);
  Double_t geomDelay=0.3e9*(1-TMath::Cos(TMath::DegToRad()*antPhi))/299792458.;
  return antDelay+geomDelay;
}

Double_t singleDelayWithGeomMod(Double_t *x, Double_t *par) {
  Double_t antPhi=x[0];
  Double_t antDelay=singleDelayFuncMod(x,par);
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

Double_t groupDelayFuncMod(Double_t *x, Double_t *par) {
  Double_t leftPhi=x[0]+11.25;
  Double_t rightPhi=x[0]-11.25;
  
  //leftPhi*=-1;
  //  rightPhi*=-1;
  Double_t leftDelay=singleDelayFuncMod(&leftPhi,par);
  Double_t rightDelay=singleDelayFuncMod(&rightPhi,par);

  return (leftDelay-rightDelay);
}


void groupDelayAnita3(TH2D* histDtPhi, TH2D* histDtPhiAnt[48])
{

   TF1 *fitDelay = new TF1("fitDelay",groupDelayFuncMod,-90,90,13);
   fitDelay->SetParLimits(0,-1,-1);
   fitDelay->SetParLimits(1,-1,-1);
   fitDelay->SetParLimits(3,-1,-1);
   fitDelay->SetParLimits(5,-1,-1);
   fitDelay->SetParLimits(7,-1,-1);
   fitDelay->SetParLimits(9,-1,-1);
   fitDelay->SetParLimits(11,-1,-1);
   // fitDelay->SetParLimits(13,-1,-1);
   // fitDelay->SetParLimits(15,-1,-1);
   // fitDelay->SetParLimits(17,-1,-1);
   // fitDelay->SetParLimits(19,-1,-1);

   TCanvas *canData2 = new TCanvas("canData2","canData2",600,600); 
   canData2->Divide(1,2);
   canData2->cd(1);
   histDtPhi->Draw("colz");
   histDtPhi->GetXaxis()->SetTitle("Phi Wave - Middle Phi (Degrees)");
   histDtPhi->GetYaxis()->SetTitle("#Delta#Deltat");
   histDtPhi->GetYaxis()->SetRangeUser(-0.2,0.2);
   canData2->cd(2);
   TProfile *histDtProf = histDtPhi->ProfileX();
   histDtProf->Draw();
   histDtProf->Fit("fitDelay", "r");
   histDtProf->GetXaxis()->SetTitle("Phi Wave - Middle Phi (Degrees)");
   histDtProf->GetYaxis()->SetTitle("#Delta#Deltat");
   histDtProf->GetYaxis()->SetRangeUser(-0.1,0.1);

   
  TCanvas *canDelay = new TCanvas("canDelay","canDelay",600,600);  
  canDelay->Divide(1,2);
  {
    canDelay->cd(1);
    TF1 *antDelay = new TF1("antDelay",singleDelayFuncMod,-100,100,5);
    antDelay->SetParameters(fitDelay->GetParameters());
    antDelay->Draw();    
    antDelay->GetXaxis()->SetTitle("Phi from Boresight");
    antDelay->GetYaxis()->SetTitle("Delay (ns)");
  }
  {
    gStyle->SetFitFormat("6.3e");
    canDelay->cd(2);
    TF1 *antDelay = new TF1("antDelayGeom",singleDelayWithGeomMod,-100,100,5);
    antDelay->SetParameters(fitDelay->GetParameters());
    antDelay->Draw();    
    antDelay->GetXaxis()->SetTitle("Phi from Boresight");
    antDelay->GetYaxis()->SetTitle("Delay + Geom (ns)");
  }
  


  TCanvas *canAnt = new TCanvas("canAnt","canAnt",800,800);
  canAnt->Divide(6,8,0,0);
  for(int ant=0;ant<MAX_ANTENNAS;ant++) {
    canAnt->cd(ant+1);
    histDtPhiAnt[ant]->Draw("colz");
    histDtPhiAnt[ant]->GetXaxis()->SetTitle("Middle Phi (Degrees)");
    histDtPhiAnt[ant]->GetYaxis()->SetTitle("#Delta#Deltat");
    histDtPhiAnt[ant]->GetYaxis()->SetRangeUser(-0.2,0.2);
  }
 

  TCanvas *canProf = new TCanvas("canProf","canProf",800,800);
  canProf->Divide(6,8,0,0);
  for(int ant=0;ant<MAX_ANTENNAS;ant++) {
    canProf->cd(ant+1);
    cout << "ANTENNA :  " << ant << endl;

    TProfile *histDtProf = histDtPhiAnt[ant]->ProfileX();
    histDtProf->Draw();
    histDtProf->Fit("fitDelay", "r");
    histDtProf->GetXaxis()->SetTitle("Middle Phi (Degrees)");
    histDtProf->GetYaxis()->SetTitle("#Delta#Deltat");
    histDtProf->GetYaxis()->SetRangeUser(-0.1,0.1);
  }
  
}
