//////////////////////////////////////////////////////////////////////////////
/////  PrettyAnitaEvent.cxx        ANITA event reading class                  /////////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for plotting event stuff like waveforms and correlations/////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//ANITA Includes
#include "PrettyAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "FFTtools.h"
#include "FFTWComplex.h"


//ROOT Includes
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TVirtualFFT.h"
#include "TMinuit.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"


ClassImp(PrettyAnitaEvent);


//My incredibly dodgy approach to fitting with MINUIT that I'm going to call from within itself
// this is horrible and dangerous and should never bve done, but hey ho there we go.
void CorSumFCN(Int_t& npar, Double_t*gin,
                       Double_t&f, Double_t* par, Int_t flag)
{
   //par[0] is phiWave
   //par[1] is thetaWave
   //par[2] is numAnts (11 or 19)

  CorrelationSummary* myDodgyCorSumPtr = dynamic_cast<CorrelationSummary*>(gMinuit->GetObjectFit());
  f=myDodgyCorSumPtr->getChiSquared(par[0],par[1],11);
  //  std::cout << par[0] << "\t" << par[1] << "\t" << f << std::endl;
}



PrettyAnitaEvent::PrettyAnitaEvent(RawAnitaEvent *eventPtr,WaveCalType::WaveCalType_t calType, PrettyAnitaHk *theHk):UsefulAnitaEvent(eventPtr,calType,theHk) { }


TCanvas *PrettyAnitaEvent::getSixWaveformCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canSixWave");
    if(!can) {
      can = new TCanvas("canSixWave","canSixWave",800,600);
    }
  }
  
  can->Clear();
  can->Divide(3,2);
  
  int topAnts[3];
  int bottomAnts[3];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);

  for(int i=0;i<3;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;
     sprintf(graphTitle,"Ant: %d",topAnts[i]);
    can->cd(i+1);
    int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
    TGraph *grTop = getGraph(ciTop);
    grTop->SetLineColor(getPrettyColour(i));
    grTop->SetTitle(graphTitle);
    grTop->Draw("al");


    can->cd(i+4);
    sprintf(graphTitle,"Ant: %d",bottomAnts[i]);
    int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);
    TGraph *grBottom = getGraph(ciBottom);
    grBottom->SetLineColor(getPrettyColour(i));
    grBottom->SetTitle(graphTitle);
    grBottom->Draw("al");
    
  }
  return can;
}


TCanvas *PrettyAnitaEvent::getTenWaveformCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canTenWave");
    if(!can) {
      can = new TCanvas("canTenWave","canTenWave",1000,600);
    }
  }
  
  can->Clear();
  can->Divide(5,2);
  
  int topAnts[3];
  int bottomAnts[3];
  int nextFourAnts[4];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);
  fillNextFourAntArrays(ant,nextFourAnts);

  int newTopAnts[5]={nextFourAnts[0],topAnts[0],topAnts[1],topAnts[2],nextFourAnts[1]};
  int newBottomAnts[5]={nextFourAnts[2],bottomAnts[0],bottomAnts[1],bottomAnts[2],nextFourAnts[3]};
  

  for(int i=0;i<5;i++) {
    //    std::cout << i << "\t" << newTopAnts[i] << "\t" << bottomAnts[i] << std::endl;
     sprintf(graphTitle,"Ant: %d",newTopAnts[i]);
    can->cd(i+1);
    int ciTop=AnitaGeomTool::getChanIndexFromAntPol(newTopAnts[i],pol);
    TGraph *grTop = getGraph(ciTop);
    grTop->SetLineColor(getPrettyColour(i));
    grTop->SetTitle(graphTitle);
    grTop->Draw("al");


    can->cd(i+6);
    sprintf(graphTitle,"Ant: %d",newBottomAnts[i]);
    int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(newBottomAnts[i],pol);
    TGraph *grBottom = getGraph(ciBottom);
    grBottom->SetLineColor(getPrettyColour(i));
    grBottom->SetTitle(graphTitle);
    grBottom->Draw("al");
    
  }
  return can;
}




TCanvas *PrettyAnitaEvent::getSixFFTPowerCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canSixFFTPower");
    if(!can) {
      can = new TCanvas("canSixFFTPower","canSixFFTPower",800,600);
    }
  }
  
  can->Clear();
  can->Divide(3,2);
  
  int topAnts[3];
  int bottomAnts[3];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);

  for(int i=0;i<3;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;
     sprintf(graphTitle,"Ant: %d",topAnts[i]);
    can->cd(i+1);
    int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
    TGraph *grTop = getFFTMagnitude(ciTop);
    grTop->SetLineColor(getPrettyColour(i));
    grTop->SetTitle(graphTitle);
    grTop->Draw("al");


    can->cd(i+4);
    sprintf(graphTitle,"Ant: %d",bottomAnts[i]);
    int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);
    TGraph *grBottom = getFFTMagnitude(ciBottom);
    grBottom->SetLineColor(getPrettyColour(i));
    grBottom->SetTitle(graphTitle);
    grBottom->Draw("al");
    
  }
  return can;
}


TCanvas *PrettyAnitaEvent::getSixPowerEnvelopeCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canSixPower");
    if(!can) {
      can = new TCanvas("canSixPower","canSixPower",800,600);
    }
  }
  
  can->Clear();
  can->Divide(3,2);
  
  int topAnts[3];
  int bottomAnts[3];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);

  for(int i=0;i<3;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;
     sprintf(graphTitle,"Ant: %d",topAnts[i]);
    can->cd(i+1);
    int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
    TGraph *grTop = getSimplePowerEnvelopeGraph(ciTop);
    grTop->SetLineColor(getPrettyColour(i));
    grTop->SetTitle(graphTitle);
    grTop->Draw("al");


    can->cd(i+4);
    sprintf(graphTitle,"Ant: %d",bottomAnts[i]);
    int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);
    TGraph *grBottom = getSimplePowerEnvelopeGraph(ciBottom);
    grBottom->SetLineColor(getPrettyColour(i));
    grBottom->SetTitle(graphTitle);
    grBottom->Draw("al");
    
  }
  return can;
}


TCanvas *PrettyAnitaEvent::getSixCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canSixCorr");
    if(!can) {
      can = new TCanvas("canSixCorr","canSixCorr",800,600);
    }
  }
  
  can->Clear();
  can->Divide(3,2);
  
  int topAnts[3];
  int bottomAnts[3];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);

  for(int i=0;i<3;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;
     sprintf(graphTitle,"Ants: %d - %d",topAnts[i],bottomAnts[i]);
    can->cd(i+1);
    int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
    int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);

    TGraph *grTop = getCorrelation(ciTop,ciBottom);
    grTop->SetLineColor(getPrettyColour(i));
    grTop->SetTitle(graphTitle);
    grTop->Draw("al");


    can->cd(i+4);
    sprintf(graphTitle,"Ants: %d - %d",bottomAnts[i],topAnts[i]);
    TGraph *grBottom = getCorrelation(ciBottom,ciTop);
    grBottom->SetLineColor(getPrettyColour(i));
    grBottom->SetTitle(graphTitle);
    grBottom->Draw("al");
    
  }
  return can;
}



TCanvas *PrettyAnitaEvent::getTenCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canTenCorr");
    if(!can) {
      can = new TCanvas("canTenCorr","canTenCorr",1000,600);
    }
  }
  
  can->Clear();
  can->Divide(5,2);
  
  int topAnts1[3];
  int bottomAnts1[3];
  int nextFourAnts[4];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts1,bottomAnts1);
  fillNextFourAntArrays(ant,nextFourAnts);

  int topAnts[5]={nextFourAnts[0],topAnts1[0],topAnts1[1],topAnts1[2],nextFourAnts[1]};
  int bottomAnts[5]={nextFourAnts[2],bottomAnts1[0],bottomAnts1[1],bottomAnts1[2],nextFourAnts[3]};
  

  for(int i=0;i<5;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;
     sprintf(graphTitle,"Ants: %d - %d",topAnts[i],bottomAnts[i]);
    can->cd(i+1);
    int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
    int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);

    TGraph *grTop = getCorrelation(ciTop,ciBottom);
    grTop->SetLineColor(getPrettyColour(i));
    grTop->SetTitle(graphTitle);
    grTop->Draw("al");


    can->cd(i+4);
    sprintf(graphTitle,"Ants: %d - %d",bottomAnts[i],topAnts[i]);
    TGraph *grBottom = getCorrelation(ciBottom,ciTop);
    grBottom->SetLineColor(getPrettyColour(i));
    grBottom->SetTitle(graphTitle);
    grBottom->Draw("al");
    
  }
  return can;
}


TCanvas *PrettyAnitaEvent::getElevenCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canElevenCorr");
    if(!can) {
      can = new TCanvas("canElevenCorr","canElevenCorr",1000,400);
    }
  }
  
  can->Clear();
  can->Divide(4,3);
  
  int topAnts[3];
  int bottomAnts[3];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);


  int ciLeftTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[0],pol);
  int ciMidTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[1],pol);
  int ciRightTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[2],pol);

  int ciLeftBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[0],pol);
  int ciMidBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[1],pol);
  int ciRightBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[2],pol);
  //First row
  //  TPad *paddy = (TPad*)can->cd(1);
  //  paddy->Divide(4,1);
  {
     can->cd(1);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[1],bottomAnts[0]);
     TGraph *grL=getCorrelation(ciMidTop,ciLeftBottom);
     grL->SetLineColor(getPrettyColour(1));
     grL->SetTitle(graphTitle);
     grL->Draw("al");


     can->cd(2);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[1],topAnts[0]);
     TGraph *grML=getCorrelation(ciMidTop,ciLeftTop);
     grML->SetLineColor(getPrettyColour(2));
     grML->SetTitle(graphTitle);
     grML->Draw("al");

     can->cd(3);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[1],topAnts[2]);
     TGraph *grMR=getCorrelation(ciMidTop,ciRightTop);
     grMR->SetLineColor(getPrettyColour(3));
     grMR->SetTitle(graphTitle);
     grMR->Draw("al");

     can->cd(4);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[1],bottomAnts[2]);
     TGraph *grR=getCorrelation(ciMidTop,ciRightBottom);
     grR->SetLineColor(getPrettyColour(4));
     grR->SetTitle(graphTitle);
     grR->Draw("al");

  }

  
  //Second row
  //  TPad *paddy2 = (TPad*)can->cd(1);
  //  paddy2->Divide(3,1);

  for(int i=0;i<3;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;

     can->cd(i+5);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[i],bottomAnts[i]);
     int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
     int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);
     
     TGraph *grTop = getCorrelation(ciTop,ciBottom);
     grTop->SetLineColor(getPrettyColour(i));
     grTop->SetTitle(graphTitle);
     grTop->Draw("al");
    
  }

  //Third row
  //  TPad *paddy3 = (TPad*) can->cd(1);
  //  paddy3->Divide(4,1);
  {
     can->cd(9);
     sprintf(graphTitle,"Ants: %d - %d",bottomAnts[1],topAnts[0]);
     TGraph *grL=getCorrelation(ciMidBottom,ciLeftTop);
     grL->SetLineColor(getPrettyColour(1));
     grL->SetTitle(graphTitle);
     grL->Draw("al");


     can->cd(10);
     sprintf(graphTitle,"Ants: %d - %d",bottomAnts[1],bottomAnts[0]);
     TGraph *grML=getCorrelation(ciMidBottom,ciLeftBottom);
     grML->SetLineColor(getPrettyColour(2));
     grML->SetTitle(graphTitle);
     grML->Draw("al");

     can->cd(11);
     sprintf(graphTitle,"Ants: %d - %d",bottomAnts[1],bottomAnts[2]);
     TGraph *grMR=getCorrelation(ciMidBottom,ciRightBottom);
     grMR->SetLineColor(getPrettyColour(3));
     grMR->SetTitle(graphTitle);
     grMR->Draw("al");

     can->cd(12);
     sprintf(graphTitle,"Ants: %d - %d",bottomAnts[1],topAnts[2]);
     TGraph *grR=getCorrelation(ciMidBottom,ciRightTop);
     grR->SetLineColor(getPrettyColour(4));
     grR->SetTitle(graphTitle);
     grR->Draw("al");

  }

  return can;
}



TCanvas *PrettyAnitaEvent::getElevenInterpolationCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, Double_t deltaT, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canElevenIntCorr");
    if(!can) {
      can = new TCanvas("canElevenIntCorr","canElevenIntCorr",1000,400);
    }
  }
  
  can->Clear();
  can->Divide(4,3);
  
  int topAnts[3];
  int bottomAnts[3];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);


  int ciLeftTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[0],pol);
  int ciMidTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[1],pol);
  int ciRightTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[2],pol);

  int ciLeftBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[0],pol);
  int ciMidBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[1],pol);
  int ciRightBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[2],pol);
  //First row
  //  TPad *paddy = (TPad*)can->cd(1);
  //  paddy->Divide(4,1);
  {
     can->cd(1);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[1],bottomAnts[0]);
     TGraph *grL=getCorrelationInterpolated(ciMidTop,ciLeftBottom,deltaT);
     grL->SetLineColor(getPrettyColour(1));
     grL->SetTitle(graphTitle);
     grL->Draw("al");


     can->cd(2);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[1],topAnts[0]);
     TGraph *grML=getCorrelationInterpolated(ciMidTop,ciLeftTop,deltaT);
     grML->SetLineColor(getPrettyColour(2));
     grML->SetTitle(graphTitle);
     grML->Draw("al");

     can->cd(3);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[1],topAnts[2]);
     TGraph *grMR=getCorrelationInterpolated(ciMidTop,ciRightTop,deltaT);
     grMR->SetLineColor(getPrettyColour(3));
     grMR->SetTitle(graphTitle);
     grMR->Draw("al");

     can->cd(4);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[1],bottomAnts[2]);
     TGraph *grR=getCorrelationInterpolated(ciMidTop,ciRightBottom,deltaT);
     grR->SetLineColor(getPrettyColour(4));
     grR->SetTitle(graphTitle);
     grR->Draw("al");

  }

  
  //Second row
  //  TPad *paddy2 = (TPad*)can->cd(1);
  //  paddy2->Divide(3,1);

  for(int i=0;i<3;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;

     can->cd(i+5);
     sprintf(graphTitle,"Ants: %d - %d",topAnts[i],bottomAnts[i]);
     int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
     int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);
     
     TGraph *grTop = getCorrelationInterpolated(ciTop,ciBottom,deltaT);
     grTop->SetLineColor(getPrettyColour(i));
     grTop->SetTitle(graphTitle);
     grTop->Draw("al");
    
  }

  //Third row
  //  TPad *paddy3 = (TPad*) can->cd(1);
  //  paddy3->Divide(4,1);
  {
     can->cd(9);
     sprintf(graphTitle,"Ants: %d - %d",bottomAnts[1],topAnts[0]);
     TGraph *grL=getCorrelationInterpolated(ciMidBottom,ciLeftTop,deltaT);
     grL->SetLineColor(getPrettyColour(1));
     grL->SetTitle(graphTitle);
     grL->Draw("al");


     can->cd(10);
     sprintf(graphTitle,"Ants: %d - %d",bottomAnts[1],bottomAnts[0]);
     TGraph *grML=getCorrelationInterpolated(ciMidBottom,ciLeftBottom,deltaT);
     grML->SetLineColor(getPrettyColour(2));
     grML->SetTitle(graphTitle);
     grML->Draw("al");

     can->cd(11);
     sprintf(graphTitle,"Ants: %d - %d",bottomAnts[1],bottomAnts[2]);
     TGraph *grMR=getCorrelationInterpolated(ciMidBottom,ciRightBottom,deltaT);
     grMR->SetLineColor(getPrettyColour(3));
     grMR->SetTitle(graphTitle);
     grMR->Draw("al");

     can->cd(12);
     sprintf(graphTitle,"Ants: %d - %d",bottomAnts[1],topAnts[2]);
     TGraph *grR=getCorrelationInterpolated(ciMidBottom,ciRightTop,deltaT);
     grR->SetLineColor(getPrettyColour(4));
     grR->SetTitle(graphTitle);
     grR->Draw("al");

  }

  return can;
}


TCanvas *PrettyAnitaEvent::getSixInterpolatedCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, Double_t deltaT, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canSixIntCorr");
    if(!can) {
      can = new TCanvas("canSixIntCorr","canSixIntCorr",800,600);
    }
  }
  
  can->Clear();
  can->Divide(3,2);
  
  int topAnts[3];
  int bottomAnts[3];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);

  for(int i=0;i<3;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;
     sprintf(graphTitle,"Ants: %d - %d",topAnts[i],bottomAnts[i]);
    can->cd(i+1);
    int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
    int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);

    TGraph *grTop = getCorrelationInterpolated(ciTop,ciBottom,deltaT);
    grTop->SetLineColor(getPrettyColour(i));
    grTop->SetTitle(graphTitle);
    grTop->Draw("al");


    can->cd(i+4);
    sprintf(graphTitle,"Ants: %d - %d",bottomAnts[i],topAnts[i]);
    TGraph *grBottom = getCorrelationInterpolated(ciBottom,ciTop,deltaT);
    grBottom->SetLineColor(getPrettyColour(i));
    grBottom->SetTitle(graphTitle);
    grBottom->Draw("al");
    
  }
  return can;
}


TCanvas *PrettyAnitaEvent::getSixInterpolatedCanvas(int ant, AnitaPol::AnitaPol_t pol, Double_t deltaT, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>31) {
    std::cerr << "Antenna " << ant << " is not in range 0-31" << std::endl;
    return NULL;
  }
  setStyleSixCanvas();
  
  //Make Canvas
  if(!can) {
    can = (TCanvas*) gROOT->FindObject("canSixInter");
    if(!can) {
      can = new TCanvas("canSixInter","canSixInter",800,600);
    }
  }
  
  can->Clear();
  can->Divide(3,2);
  
  int topAnts[3];
  int bottomAnts[3];
  //  std::cout << ant << "\t" << AnitaGeomTool::getAzimuthPartner(ant) << std::endl;
  fillSixAntArrays(ant,topAnts,bottomAnts);

  for(int i=0;i<3;i++) {
    //    std::cout << i << "\t" << topAnts[i] << "\t" << bottomAnts[i] << std::endl;
     sprintf(graphTitle,"Ant: %d",topAnts[i]);
     can->cd(i+1);
     int ciTop=AnitaGeomTool::getChanIndexFromAntPol(topAnts[i],pol);
     TGraph *grTop = getInterpolatedGraph(ciTop,deltaT);
     grTop->SetLineColor(getPrettyColour(i));
     grTop->SetTitle(graphTitle);
     grTop->Draw("al");
     
     
     can->cd(i+4);
     sprintf(graphTitle,"Ant: %d",bottomAnts[i]);
     int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(bottomAnts[i],pol);
     TGraph *grBottom = getInterpolatedGraph(ciBottom,deltaT);
     grBottom->SetLineColor(getPrettyColour(i));
     grBottom->SetTitle(graphTitle);
     grBottom->Draw("al");
     
  }
  return can; 
}


TGraph *PrettyAnitaEvent::getSimplePowerEnvelopeGraph(int chanIndex) {
  //Need to loop over each channel and form V^2 array and then loop
  // over that array and pick only the channels that are larger than (or equal to their neighbours.
  Double_t v2[NUM_SAMP];
  Double_t vEnvelope[NUM_SAMP];
  Double_t tEnvelope[NUM_SAMP];
  Int_t numPoints=0;

  for(int i=0;i<fNumPoints[chanIndex];i++) {
    v2[i]=fVolts[chanIndex][i]*fVolts[chanIndex][i];
    if(i==1) {
      if(v2[0]>v2[i]) {
	vEnvelope[numPoints]=v2[0];
	tEnvelope[numPoints]=fTimes[chanIndex][0];
	numPoints++;
      }
    }
    else if(i==fNumPoints[chanIndex]-1 && v2[i]>v2[i-1]) {
      vEnvelope[numPoints]=v2[i];
      tEnvelope[numPoints]=fTimes[chanIndex][i];
      numPoints++;
    }
    else if(v2[i-1]>v2[i-2] && v2[i-1]>v2[i]) {
      vEnvelope[numPoints]=v2[i-1];
      tEnvelope[numPoints]=fTimes[chanIndex][i-1];
      numPoints++;
    }
  }				     
  TGraph *grPower = new TGraph(numPoints,tEnvelope,vEnvelope);
  return grPower;
      
}


TGraph *PrettyAnitaEvent::getInterpolatedGraph(int chanIndex, double deltaT) {
   //Will use the ROOT::Math::Interpolator function to do this.
   std::vector<double> tVec;
   std::vector<double> vVec;

   for (int samp=0;samp<fNumPoints[chanIndex];samp++) {
     tVec.push_back(fTimes[chanIndex][samp]);
     vVec.push_back(fVolts[chanIndex][samp]);
   }
   if(tVec.size()==0) {
     std::cout << "\n" << eventNumber << "\t" << chanIndex << "\t" << tVec.size() << "\t" << vVec.size() << "\t" << fNumPoints[chanIndex] << "\n";   
   }

   ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);
   Double_t startTime=fTimes[chanIndex][0];
   Double_t lastTime=fTimes[chanIndex][fNumPoints[chanIndex]-1];


   Double_t newTimes[8192];
   Double_t newVolts[8192];
   Int_t numPoints=0;
   for(Double_t time=startTime;time<=lastTime;time+=deltaT) {
      newTimes[numPoints]=time;
      newVolts[numPoints]=chanInterp.Eval(time);
      //      std::cout << numPoints << "\t" << newTimes[numPoints]
      //		<< "\t" << newVolts[numPoints] << std::endl;
	       
      numPoints++;
   }

   TGraph *grInt = new TGraph(numPoints,newTimes,newVolts);
   return grInt;
}

TGraph *PrettyAnitaEvent::getFFTMagnitude(int chanIndex)
{
   TGraph *grWave=getGraph(chanIndex);
   TGraph *grFFT=getFFTMagnitude(grWave);
   delete grWave;
   return grFFT;
}


TGraph *PrettyAnitaEvent::getFFTMagnitude(TGraph *grIn)
{
   TGraph *grFFT=FFTtools::makePowerSpectrum(grIn);
   return grFFT;
}

TGraph *PrettyAnitaEvent::getCorrelation(TGraph *gr1, TGraph *gr2)
{
   return FFTtools::getCorrelationGraph(gr1,gr2);
}

TGraph *PrettyAnitaEvent::getCorrelation(int chanIndex1, int chanIndex2) 
{
   TGraph *gr1 =getGraph(chanIndex1);
   TGraph *gr2 = getGraph(chanIndex2);
   TGraph *grCor = getCorrelation(gr1,gr2);
   Double_t x1,y1,x2,y2;
   gr1->GetPoint(0,x1,y1);
   //   std::cout << 1 << "\t" << chanIndex1 << "\t" << x << "\t" << y << "\n";
   gr2->GetPoint(0,x2,y2);
   fWaveOffset=x1-x2;
   
   //   std::cout << 2 << "\t" << chanIndex2 << "\t" << x << "\t" << y << "\n";   
   //   grCor->GetPoint(0,x1,y1);
   gr1->GetPoint(1,x2,y2);
   fDeltaT=x2-x1;
   
   //   Double_t x,y;
   //   gr1->GetPoint(0,x,y);
   //   std::cout << 1 << "\t" << chanIndex1 << "\t" << x << "\t" << y << "\n";
   //   gr2->GetPoint(0,x,y);
   //   std::cout << 2 << "\t" << chanIndex2 << "\t" << x << "\t" << y << "\n";
   //   grCor->GetPoint(0,x,y);
   //   std::cout << "Cor\t" << chanIndex2 << "\t" << x << "\t" << y << "\n";


   delete gr1;
   delete gr2;
   return grCor;     
}

TGraph *PrettyAnitaEvent::getCorrelationInterpolated(int chanIndex1, int chanIndex2, Double_t deltaT) 
{
  if(chanIndex1 <0 || chanIndex1>80)
    std::cerr << "Invalid channel index:\t" << chanIndex1 << "\n";
  if(chanIndex2 <0 || chanIndex2>80)
    std::cerr << "Invalid channel index:\t" << chanIndex2 << "\n";
   TGraph *gr1 =getInterpolatedGraph(chanIndex1,deltaT);
   TGraph *gr2 = getInterpolatedGraph(chanIndex2,deltaT);
   TGraph *grCor = getCorrelation(gr1,gr2);
   Double_t x1,y1,x2,y2;
   gr1->GetPoint(0,x1,y1);
   //   std::cout << 1 << "\t" << chanIndex1 << "\t" << x << "\t" << y << "\n";
   gr2->GetPoint(0,x2,y2);
   fWaveOffset=x1-x2;
   
   //   std::cout << 2 << "\t" << chanIndex2 << "\t" << x << "\t" << y << "\n";   
   //   grCor->GetPoint(0,x1,y1);
   gr1->GetPoint(1,x2,y2);
   fDeltaT=x2-x1;
   
   //   std::cout << "CorInt\t" << chanIndex2 << "\t" << x << "\t" << y << "\n";
   delete gr1;
   delete gr2;
   return grCor;
}




void PrettyAnitaEvent::setStyleSixCanvas()
{
   gStyle->SetPadTopMargin(0.05);
   gStyle->SetPadBottomMargin(0.05);
   gStyle->SetPadLeftMargin(0.05);
   gStyle->SetPadRightMargin(0.05);
}

void PrettyAnitaEvent::fillSixAntArrays(int ant, int topAnts[3], int bottomAnts[3])
{
   //   std::cerr << "fillSixAntArrays( " << ant << ")\n";
  int top=-1,bottom=-1;
  int leftTop=-1, rightTop=-1;
  int leftBottom=-1, rightBottom=-1;

  if(ant<16) {
    top=ant;
    bottom=AnitaGeomTool::getAzimuthPartner(top);
  }
  else {
    bottom=ant;
    top=AnitaGeomTool::getAzimuthPartner(bottom);
  }
  //  std::cout << top << "\t" << bottom << std::endl;
  AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
  AnitaGeomTool::getThetaPartners(bottom,leftBottom,rightBottom);
  topAnts[0]=leftTop;
  topAnts[1]=top;
  topAnts[2]=rightTop;
  bottomAnts[0]=leftBottom;
  bottomAnts[1]=bottom;
  bottomAnts[2]=rightBottom;

}

void PrettyAnitaEvent::fillNextFourAntArrays(int ant, int nextFourAnts[4])
{

  int top=-1,bottom=-1;
  int leftTop=-1, rightTop=-1;
  int leftLeftTop=-1, rightRightTop=-1;
  int leftBottom=-1, rightBottom=-1;
  int leftLeftBottom=-1, rightRightBottom=-1;

  if(ant<16) {
    top=ant;
    bottom=AnitaGeomTool::getAzimuthPartner(top);
  }
  else {
    bottom=ant;
    top=AnitaGeomTool::getAzimuthPartner(bottom);
  }
  int crap;
  //  std::cout << top << "\t" << bottom << std::endl;
  AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
  AnitaGeomTool::getThetaPartners(bottom,leftBottom,rightBottom);
  AnitaGeomTool::getThetaPartners(leftTop,leftLeftTop,crap);
  AnitaGeomTool::getThetaPartners(rightTop,crap,rightRightTop);
  AnitaGeomTool::getThetaPartners(leftBottom,leftLeftBottom,crap);
  AnitaGeomTool::getThetaPartners(rightBottom,crap,rightRightBottom);
  nextFourAnts[0]=leftLeftTop;
  nextFourAnts[1]=rightRightTop;
  nextFourAnts[2]=leftLeftBottom;
  nextFourAnts[3]=rightRightBottom;

}

int PrettyAnitaEvent::getPrettyColour(int index)
{
    if(index>10) return index;
    Int_t niceColours[11]={50,38,30,9,8,44,24,12,40,20,41};
    return niceColours[index];
}


//Putative Analysis methods
int PrettyAnitaEvent::getMaxAntenna(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
   return getMaxAntennaCorrelation(pol, peakPtr);
}   


int PrettyAnitaEvent::getMaxAntennaVSquared(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
   //Returns the antenna with the maximum power
   //Could consider changng this to make things better
   double maxVal=0;
   int maxAnt=0;
   for(int ant=0;ant<32;ant++) {
      int chanIndex=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
      for(int samp=0;samp<fNumPoints[chanIndex];samp++) {
	 double vSquared=fVolts[chanIndex][samp]*fVolts[chanIndex][samp];
	 if(vSquared>maxVal) {
	    maxVal=vSquared;
	    maxAnt=ant;
	 }
      }
   }
   if(peakPtr) *peakPtr=maxVal;
   return maxAnt;
}

int PrettyAnitaEvent::getMaxAntennaCorrelation(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
   //Returns the antenna with the lagest peak/rms in the correlation with its azimuth partner antenna
   double maxVal=0;
   int maxAnt=0;
   for(int ant=0;ant<16;ant++) {
      //Loop over the top antennas
      int otherAnt=AnitaGeomTool::getAzimuthPartner(ant);
      int ciTop=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
      int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(otherAnt,pol);

      TGraph *grCor = getCorrelation(ciTop,ciBottom);      
      Double_t *y = grCor->GetY();
      Double_t peak=TMath::MaxElement(grCor->GetN(),y);
      Double_t rms=TMath::RMS(grCor->GetN(),y);
      //      FFTtools::getPeakRmsSqVal(grCor,peak,rms);
      if((peak/rms)>maxVal) {
	 maxVal=peak/rms;
	 maxAnt=ant;
	 Double_t maxTop=TMath::MaxElement(fNumPoints[ciTop],fVolts[ciTop]);
	 Double_t maxBottom=TMath::MaxElement(fNumPoints[ciBottom],fVolts[ciBottom]);
	 if(maxBottom>maxTop)
	   maxAnt=otherAnt;
      }
      delete grCor;
	 
   }
   if(peakPtr) *peakPtr=maxVal;
   return maxAnt;
}


CorrelationSummary *PrettyAnitaEvent::getCorrelationSummary(Int_t centreAnt,AnitaPol::AnitaPol_t pol, Double_t deltaT)
{
  //Gets the 11 correlations and then takes the max, rms and neighbouring maxima
  if(centreAnt<0)
    centreAnt=getMaxAntenna(pol);
  int sixAnts[6];
  fillSixAntArrays(centreAnt,sixAnts,&(sixAnts[3]));
  int nextFourAnts[4];
  fillNextFourAntArrays(centreAnt,nextFourAnts);


     
   
   CorrelationSummary *theSum = new CorrelationSummary(eventNumber, centreAnt, sixAnts,deltaT);
   for(int i=0;i<4;i++)
      theSum->nextFourAnts[i]=nextFourAnts[i];


   //Now need to make correlation index pairs
   //Top-Bottom first
   theSum->firstAnt[0]=sixAnts[0];
   theSum->secondAnt[0]=sixAnts[3];
   theSum->firstAnt[1]=sixAnts[1];
   theSum->secondAnt[1]=sixAnts[4];
   theSum->firstAnt[2]=sixAnts[2];
   theSum->secondAnt[2]=sixAnts[5];
   //Now Left-Right
   theSum->firstAnt[3]=sixAnts[0];
   theSum->secondAnt[3]=sixAnts[1];
   theSum->firstAnt[4]=sixAnts[1];
   theSum->secondAnt[4]=sixAnts[2];
   theSum->firstAnt[5]=sixAnts[3];
   theSum->secondAnt[5]=sixAnts[4];
   theSum->firstAnt[6]=sixAnts[4];
   theSum->secondAnt[6]=sixAnts[5];
   //Now Diagonal
   theSum->firstAnt[7]=sixAnts[0];
   theSum->secondAnt[7]=sixAnts[4];
   theSum->firstAnt[8]=sixAnts[2];
   theSum->secondAnt[8]=sixAnts[4];
   theSum->firstAnt[9]=sixAnts[1];
   theSum->secondAnt[9]=sixAnts[3];
   theSum->firstAnt[10]=sixAnts[1];
   theSum->secondAnt[10]=sixAnts[5];


   //Now Leftmost - centre top
   theSum->firstAnt[11]=nextFourAnts[0];
   theSum->secondAnt[11]=sixAnts[1];
   //Now centre - right most top
   theSum->firstAnt[12]=sixAnts[1];
   theSum->secondAnt[12]=nextFourAnts[1];
   //Now Leftmost - centre bottom
   theSum->firstAnt[13]=nextFourAnts[2];
   theSum->secondAnt[13]=sixAnts[4];
   //Now centre - right most bottom
   theSum->firstAnt[14]=sixAnts[4];
   theSum->secondAnt[14]=nextFourAnts[3];

   //Now Leftmost - left top
   theSum->firstAnt[15]=nextFourAnts[0];
   theSum->secondAnt[15]=sixAnts[0];
   //Now right - right most top
   theSum->firstAnt[16]=sixAnts[2];
   theSum->secondAnt[16]=nextFourAnts[1];
   //Now Leftmost - left bottom
   theSum->firstAnt[17]=nextFourAnts[2];
   theSum->secondAnt[17]=sixAnts[3];
   //Now right - right most bottom
   theSum->firstAnt[18]=sixAnts[5];
   theSum->secondAnt[18]=nextFourAnts[3];
   
   



   //Now can make correlations and find max, rms, etc.
   for(int corInd=0;corInd<19;corInd++) {
      TGraph *grCor;
      Int_t ci1=AnitaGeomTool::getChanIndexFromAntPol(theSum->firstAnt[corInd],pol);
      Int_t ci2=AnitaGeomTool::getChanIndexFromAntPol(theSum->secondAnt[corInd],pol);
      if(deltaT==0) {
	 grCor=getCorrelation(ci1,ci2);
      }
      else {
	 grCor=getCorrelationInterpolated(ci1,ci2,deltaT);
      }

      //      theSum->rmsCorVals[corInd]=grCor->GetRMS(2);

      double *theTimes = grCor->GetX();
      double *theValues = grCor->GetY();
      
      int numPoints=grCor->GetN();
      double rmsVal=TMath::RMS(numPoints,theValues);
      int maxIndex=TMath::LocMax(numPoints,theValues);;
      double maxVal=theValues[maxIndex];;
      //      Double_t maxVal,rmsVal;
      //      Int_t maxIndex;
      //      FFTtools::getPeakRmsRectified(grCor,maxVal,rmsVal,&maxIndex);

      //      FFTtools::getPeakRmsSqVal(grCor,maxVal,rmsVal,&maxIndex);
      //      for(int i=0;i<grCor->GetN();i++) {
      //	std::cout << i << "\t" << theTimes[i] << "\t" << theValues[i] << "\n";
      //	 if(theValues[i]>maxVal) {
      //	    maxVal=theValues[i];
      //	    maxIndex=i;
      //	 }
      //      }
      theSum->rmsCorVals[corInd]=rmsVal;
      theSum->maxCorVals[corInd]=theValues[maxIndex];
      theSum->maxCorTimes[corInd]=theTimes[maxIndex];

      //      std::cout << theSum->firstAnt[corInd] << "\t" << theSum->secondAnt[corInd]
// 		     << "\t" << theSum->maxCorTimes[corInd] 
// 		     << "\t" << theSum->maxCorVals[corInd] << "\t" 
// 		     << "\t" << (theSum->maxCorTimes[corInd]-fWaveOffset)/fDeltaT << "\t"
// 		     << fWaveOffset << "\t" << fDeltaT << std::endl;

      theSum->secondCorVals[corInd][0]=theSum->maxCorVals[corInd];
      theSum->secondCorTimes[corInd][0]=theSum->maxCorTimes[corInd];
      theSum->secondCorVals[corInd][1]=theSum->maxCorVals[corInd];
      theSum->secondCorTimes[corInd][1]=theSum->maxCorTimes[corInd];
      for(int i=maxIndex-1;i>=1;i--) {
	 if(i<1) break;	 
	 if(theValues[i]>=theValues[i-1] && theValues[i]>=theValues[i+1]) {
	    theSum->secondCorVals[corInd][0]=theValues[i];
	    theSum->secondCorTimes[corInd][0]=theTimes[i];
	    break;
	 }	  
      }
      for(int i=maxIndex+1;i<grCor->GetN();i++) {
	 if(i>=grCor->GetN()-1) break;	 
	 if(theValues[i]>=theValues[i-1] && theValues[i]>=theValues[i+1]) {
	    theSum->secondCorVals[corInd][1]=theValues[i];
	    theSum->secondCorTimes[corInd][1]=theTimes[i];
	    break;
	 }	  
      }   
      delete grCor;
   }

   //Will add a call to
   theSum->fillErrorsAndFit();

   //Set up MINUIT for the fit
   static int firstTime=1;
   if(firstTime) {
      gMinuit = new TMinuit(2);
      firstTime=0;
   }
   gMinuit->SetObjectFit(theSum);  
   gMinuit->SetFCN(CorSumFCN);
   double par[2]={theSum->fAntPhi[1][0],0};               // the start values
   double stepSize[2]={0.01,0.01};          // step sizes 
   double minVal[2]={0,-1*TMath::PiOver2()};            // minimum bound on parameter 
   double maxVal[2]={TMath::TwoPi(),TMath::PiOver2()};            // maximum bound on parameter
   char parName[2][20];
   sprintf(parName[0],"phiWave");
   sprintf(parName[1],"thetaWave");
   for (int i=0; i<2; i++){
      gMinuit->DefineParameter(i, parName[i], par[i], stepSize[i], minVal[i], maxVal[i]);
   }
   
   Double_t phiWave,thetaWave;
   Double_t phiWaveErr,thetaWaveErr;
   //do the fit and get the answers
   gMinuit->SetPrintLevel(-1);
   gMinuit->Migrad();       // Minuit's best minimization algorithm   
   gMinuit->GetParameter(0,phiWave,phiWaveErr);
   gMinuit->GetParameter(1,thetaWave,thetaWaveErr);

   Int_t npari,nparx,istat;
   Double_t fmin,fedm,errdef;
   gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
   //   std::cout << fmin << "\t" << fedm << "\t" << npari << "\t" << nparx 
   //	     << "\t" << istat << std::endl;
   theSum->setFitResults(phiWave,thetaWave,phiWaveErr,thetaWaveErr,fmin);



   return theSum;
}

