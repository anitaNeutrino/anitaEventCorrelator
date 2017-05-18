//////////////////////////////////////////////////////////////////////////////
/////  PrettyAnitaEvent.cxx        ANITA event reading class                  /////////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for plotting event stuff like waveforms and correlations/////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include "PrettyAnitaEvent.h"

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

//My incredibly dodgy approach to fitting with MINUIT that I'm going to call from within itself
// this is horrible and dangerous and should never bve done, but hey ho there we go.
void CorSumFCNanita3(Int_t& npar, Double_t*gin,
                       Double_t&f, Double_t* par, Int_t flag)
{
   //par[0] is phiWave
   //par[1] is thetaWave
   //par[2] is numAnts (11 or 19)

  CorrelationSummaryAnita3* myDodgyCorSumPtr = dynamic_cast<CorrelationSummaryAnita3*>(gMinuit->GetObjectFit());
  f=myDodgyCorSumPtr->getChiSquared(par[0],par[1],11);
  //  std::cout << par[0] << "\t" << par[1] << "\t" << f << std::endl;
}

PrettyAnitaEvent::PrettyAnitaEvent(CalibratedAnitaEvent *eventPtr, WaveCalType::WaveCalType_t calType):UsefulAnitaEvent(eventPtr,calType) {
fPassBandFilter=0;
fNotchFilter=0;
}

PrettyAnitaEvent::PrettyAnitaEvent(RawAnitaEvent *eventPtr,WaveCalType::WaveCalType_t calType, PrettyAnitaHk *theHk):UsefulAnitaEvent(eventPtr,calType,theHk) {
fPassBandFilter=0;
fNotchFilter=0;
}

PrettyAnitaEvent::PrettyAnitaEvent(RawAnitaEvent *eventPtr,WaveCalType::WaveCalType_t calType,RawAnitaHeader *headPtr):UsefulAnitaEvent(eventPtr,calType,headPtr) {
   fPassBandFilter=0;
   fNotchFilter=0;
}


TCanvas *PrettyAnitaEvent::getSixWaveformCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *incan) {
   char graphTitle[180];
   
  TCanvas *can=incan;
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" << (NUM_SEAVEYS-1) << std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" << (NUM_SEAVEYS-1) << std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" << (NUM_SEAVEYS-1) << std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" <<(NUM_SEAVEYS-1)<< std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" <<(NUM_SEAVEYS-1) << std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" <<(NUM_SEAVEYS-1) << std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" <<(NUM_SEAVEYS-1) << std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" << (NUM_SEAVEYS-1) <<std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" << (NUM_SEAVEYS-1) << std::endl;
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
  if(ant<0 || ant>(NUM_SEAVEYS-1)) {
    std::cerr << "Antenna " << ant << " is not in range 0-" << (NUM_SEAVEYS-1) << std::endl;
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
   TGraph *grWave = getGraph(chanIndex);
   TGraph *grPowerEnv = FFTtools::getSimplePowerEnvelopeGraph(grWave);
   delete grWave;
   return grPowerEnv;      
}


TGraph *PrettyAnitaEvent::getInterpolatedGraph(int chanIndex, double deltaT) {
   
   TGraph *grWave = getGraph(chanIndex);

   double x = 0; double x1 = 0;
   double y = 0; double y1 = 0;

   int numP =  grWave->GetN();
   
   for(int i = 1; i < numP; i++){
      grWave->GetPoint(i,x,y);
      grWave->GetPoint(i-1,x1,y1);
      
      if((x<x1)){
	 std::cout << x << "  " << x1 << "  "<< chanIndex << "  " << i << "  " << numP << "  "<< eventNumber << std::endl;
     }
      
   }

   if(fPassBandFilter || fNotchFilter) {
     TGraph *grInt = FFTtools::getInterpolatedGraph(grWave,1./2.6);
     delete grWave;     
     //Now for the filtering
     if(fPassBandFilter) {
       TGraph *grFilt = FFTtools::simplePassBandFilter(grInt,fLowPassEdge,fHighPassEdge);
       delete grInt;
       grInt=grFilt;
     }
     if(fNotchFilter>0) {
       TGraph *grFilt = FFTtools::multipleSimpleNotchFilters(grInt,fNotchFilter,fLowNotchEdge,fHighNotchEdge);
       delete grInt;
       grInt=grFilt;
     }
     TGraph *grFinal = FFTtools::getInterpolatedGraph(grInt,deltaT);
     delete grInt;
     return grFinal;
   }

   
   TGraph *grInt = FFTtools::getInterpolatedGraph(grWave,deltaT);
   delete grWave;  
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
  if(chanIndex1 <0 || chanIndex1>(NUM_DIGITZED_CHANNELS-1))
    std::cerr << "Invalid channel index:\t" << chanIndex1 << "\n";

  if(chanIndex2 <0 || chanIndex2>(NUM_DIGITZED_CHANNELS-1))
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


 //   TFile *fp = new TFile("tempCorFile.root","RECREATE");
//    grCor->SetName("grCor");
//    grCor->Write();
//    gr1->SetName("gr1");
//    gr1->Write();
//    gr2->SetName("gr2");
//    gr2->Write();
//    TGraph *fft1 = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(gr1);
//    fft1->SetName("fft1");
//    fft1->Write();
//    TGraph *fft2 = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(gr2);
//    fft2->SetName("fft2");
//    fft2->Write();

//    exit(0);

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


void PrettyAnitaEvent::fillNineAntArrays(int ant, int nineAnts[9])
{
  // Top 0-2
  // Middle 3-5
  // Bottom 6-8
  
   //   std::cerr << "fillSixAntArrays( " << ant << ")\n";
  int top=-1,middle=-1,bottom=-1;
  int leftTop=-1, rightTop=-1;
  int leftMiddle=-1, rightMiddle=-1;
  int leftBottom=-1, rightBottom=-1;

  if(ant<16) {
    top=ant;
    middle=AnitaGeomTool::getAzimuthPartner(top);
    bottom=middle+16;
  } else if (ant<32){
    middle = ant;
    top = AnitaGeomTool::getAzimuthPartner(middle);
    bottom=middle+16;
  } else {
    bottom=ant;
    middle=AnitaGeomTool::getAzimuthPartner(bottom);
    top=AnitaGeomTool::getAzimuthPartner(middle);
  } 

//   std::cout << top << "\t" << bottom << std::endl;
  AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
  AnitaGeomTool::getThetaPartners(middle,leftMiddle,rightMiddle);
  AnitaGeomTool::getThetaPartners(bottom,leftBottom,rightBottom);

  nineAnts[0] = leftTop;
  nineAnts[1] = top;
  nineAnts[2] = rightTop;
  nineAnts[3] = leftMiddle;
  nineAnts[4] = middle;
  nineAnts[5] = rightMiddle;
  nineAnts[6] = leftBottom;
  nineAnts[7] = bottom;
  nineAnts[8] = rightBottom;

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

void PrettyAnitaEvent::fillNextSixAntArrays(int ant, int nextSixAnts[6])
{

  int top=-1,middle=-1,bottom=-1;
  int leftTop=-1, rightTop=-1;
  int leftLeftTop=-1, rightRightTop=-1;
  int leftMiddle=-1, rightMiddle=-1;
  int leftLeftMiddle=-1, rightRightMiddle=-1;
  int leftBottom=-1, rightBottom=-1;
  int leftLeftBottom=-1, rightRightBottom=-1;

  if(ant<16) {
    top=ant;
    middle=AnitaGeomTool::getAzimuthPartner(top);
    bottom=middle+16;
  } else if (ant<32){
    middle = ant;
    top = AnitaGeomTool::getAzimuthPartner(middle);
    bottom=middle+16;
  } else {
    bottom=ant;
    middle=AnitaGeomTool::getAzimuthPartner(bottom);
    top=AnitaGeomTool::getAzimuthPartner(middle);
  }


  int crap;
  //  std::cout << top << "\t" << bottom << std::endl;
  AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
  AnitaGeomTool::getThetaPartners(middle,leftMiddle,rightMiddle);
  AnitaGeomTool::getThetaPartners(bottom,leftBottom,rightBottom);
  AnitaGeomTool::getThetaPartners(leftTop,leftLeftTop,crap);
  AnitaGeomTool::getThetaPartners(rightTop,crap,rightRightTop);
  AnitaGeomTool::getThetaPartners(leftMiddle,leftLeftMiddle,crap);
  AnitaGeomTool::getThetaPartners(rightMiddle,crap,rightRightMiddle);
  AnitaGeomTool::getThetaPartners(leftBottom,leftLeftBottom,crap);
  AnitaGeomTool::getThetaPartners(rightBottom,crap,rightRightBottom);
  nextSixAnts[0]=leftLeftTop;
  nextSixAnts[1]=rightRightTop;
  nextSixAnts[2]=leftLeftMiddle;
  nextSixAnts[3]=rightRightMiddle;
  nextSixAnts[4]=leftLeftBottom;
  nextSixAnts[5]=rightRightBottom;

}


void PrettyAnitaEvent::fillNadirArrays(int ant, int nadirAnts[5])
{
  //   std::cerr << "fillSixAntArrays( " << ant << ")\n";
 
  //some logic to convert from top to nadir (on or between antennas)



  int nadirAntNums[NUM_PHI]={32,-1,33,-1,34,-1,35,-1,36,-1,37,-1,38,-1,39,-1};

  int left = -1; int right = -1;
  int leftLeft = -1; int rightRight = -1;
  int antBottom = ant;  

    if(ant<16) {
      antBottom=AnitaGeomTool::getAzimuthPartner(ant);
    }
      
    int nadir =  nadirAntNums[antBottom-16];

    if(nadir == -1){
      
      leftLeft = -1; 
      rightRight = -1;
	left = nadirAntNums[antBottom-17];

	if((antBottom-15)==16){
	right = 32;
      }else{
	right = nadirAntNums[antBottom-15];
      }  

    }else{
      left = -1; 
      right = -1;
      if(nadir==32){
	leftLeft = 39;
      }else{
	leftLeft = nadir - 1;
      }

      if(nadir==39){
	rightRight = 32;
      }else{
	rightRight = nadir + 1;
      }  
    }

    //only need nadir antennas

 // AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
//   AnitaGeomTool::getThetaPartners(ant,leftBottom,rightBottom);
//   AnitaGeomTool::getThetaPartners(leftTop,leftLeftTop,crap);
//   AnitaGeomTool::getThetaPartners(rightTop,crap,rightRightTop);
//   AnitaGeomTool::getThetaPartners(leftBottom,leftLeftBottom,crap);
//   AnitaGeomTool::getThetaPartners(rightBottom,crap,rightRightBottom);

//   //  std::cout << top << "\t" << bottom << std::endl;
//   AnitaGeomTool::getThetaPartners(nadir,leftTop,rightTop);
	nadirAnts[0]=leftLeft;
	nadirAnts[1]=rightRight;
	nadirAnts[2]=left;
	nadirAnts[3]=right;
	nadirAnts[4]=nadir;
//   nadirAnts[5]=rightRightTop;
//   nadirAnts[6]=ant;
//   nadirAnts[7]=top;
  

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
  return getMaxAntennaCorrelationRollingAvg(pol, peakPtr);
  //  return getMaxAntennaCorrelation(pol, peakPtr);
  //  getMaxAntennaVSquared(pol,peakPtr);
}   


int PrettyAnitaEvent::getMaxAntennaVSquared(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
   //Returns the antenna with the maximum power
   //Could consider changng this to make things better
   double maxVal=0;
   int maxAnt=0;
   for(int ant=0;ant<32;ant++) {
      int chanIndex=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
      Double_t rmsChan=TMath::RMS(fNumPoints[chanIndex],fVolts[chanIndex]);
      for(int samp=0;samp<fNumPoints[chanIndex];samp++) {
	 double vSquared=fVolts[chanIndex][samp]*fVolts[chanIndex][samp];
	 vSquared/=rmsChan;
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

int PrettyAnitaEvent::getMaxAntennaCorrelationRollingAvg(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
   //Returns the antenna at the centre of three phi-secotrs with the largest peak/rms in the correlation with its azimuth partner antenna
   double maxVal=0;
   int maxAnt=0;
   double maxVals[16]={0};
   for(int ant=0;ant<16;ant++) {
      //Loop over the top antennas
      int otherAnt=AnitaGeomTool::getAzimuthPartner(ant);
      int ciTop=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
      int ciMiddle=AnitaGeomTool::getChanIndexFromAntPol(otherAnt,pol);
      int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(otherAnt+16,pol);

      TGraph *grCor1 = getCorrelation(ciTop,ciMiddle);      
      TGraph *grCor2 = getCorrelation(ciTop,ciBottom);      
      TGraph *grCor3 = getCorrelation(ciMiddle,ciBottom);      
      Double_t *y1 = grCor1->GetY();
      Double_t PeakRMS[3];
      PeakRMS[0]=TMath::MaxElement(grCor1->GetN(),y1)/TMath::RMS(grCor1->GetN(),y1);
      Double_t *y2 = grCor2->GetY();
      PeakRMS[1]=TMath::MaxElement(grCor2->GetN(),y2)/TMath::RMS(grCor2->GetN(),y2);

      Double_t *y3 = grCor3->GetY();
      PeakRMS[2]=TMath::MaxElement(grCor3->GetN(),y3)/TMath::RMS(grCor3->GetN(),y3);
      Double_t maxPeakRMS = TMath::MaxElement(3, PeakRMS);

      //      FFTtools::getPeakRmsSqVal(grCor,peak,rms);
      if(maxPeakRMS>maxVals[ant])
	maxVals[ant]=maxPeakRMS;

     delete grCor1;
     delete grCor2;
     delete grCor3;
   }
   
   for(int ant=0;ant<16;ant++) {
     int leftAnt=ant-1;
     if(leftAnt<0) leftAnt=15;
     int rightAnt=ant+1;
     if(rightAnt>15) rightAnt=0;


     int otherAnt=AnitaGeomTool::getAzimuthPartner(ant);
     int ciTop=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
     int ciMiddle=AnitaGeomTool::getChanIndexFromAntPol(otherAnt,pol);
     int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(otherAnt+16,pol);
     
     Double_t newVal=maxVals[leftAnt]+maxVals[ant]+maxVals[rightAnt];

     if(newVal>maxVal) {
       maxVal=newVal;
       maxAnt=ant;
       Double_t maxTop=TMath::MaxElement(fNumPoints[ciTop],fVolts[ciTop]);
       Double_t maxMiddle=TMath::MaxElement(fNumPoints[ciMiddle],fVolts[ciMiddle]);
       Double_t maxBottom=TMath::MaxElement(fNumPoints[ciBottom],fVolts[ciBottom]);
       if(maxTop>maxBottom && maxTop>maxMiddle) maxAnt = ant;
       else if (maxMiddle>maxBottom) maxAnt = otherAnt;
       else maxAnt=otherAnt+16;
     }
	 
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
  int nadirAnts[5]={0};
  fillNadirArrays(centreAnt,nadirAnts);
     
   
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
   
   
   //Now Nadir

   if(nadirAnts[4]>-1){
      //      std::cout << "Here1\t" << nadirAnts[4] << "\t" << nadirAnts[0] << "\n";
      //   -  "[19]" CN - LLN
   theSum->firstAnt[19]=nadirAnts[4];
   theSum->secondAnt[19]=nadirAnts[0];
   //   std::cout << "Here1\t" << theSum->firstAnt[19] << "\t" << theSum->secondAnt[19] << "\n";

   //   -  "[20]" CN - RRN
   theSum->firstAnt[20]=nadirAnts[4];
   theSum->secondAnt[20]=nadirAnts[1];

  //   -  "[21]" LLN - LLB
   theSum->firstAnt[21]=nadirAnts[4];
   theSum->secondAnt[21]=nextFourAnts[2];

  //  -  "[22]" LLN - LB
   theSum->firstAnt[22]=nadirAnts[0];
   theSum->secondAnt[22]=sixAnts[3];

//   -  "[23]" CN - CB
   theSum->firstAnt[23]=nadirAnts[4];
   theSum->secondAnt[23]=sixAnts[4];

 //   -  "[24]" CN - LB
   theSum->firstAnt[24]=nadirAnts[4];
   theSum->secondAnt[24]=sixAnts[3];

//   -  "[25]" CN - RB
   theSum->firstAnt[25]=nadirAnts[4];
   theSum->secondAnt[25]=sixAnts[5];

 //   -  "[26]" RRN - RB
   theSum->firstAnt[26]=nadirAnts[1];
   theSum->secondAnt[26]=sixAnts[5];

 //   -  "[27]" RRN - RRB
   theSum->firstAnt[27]=nadirAnts[1];
   theSum->secondAnt[27]=nextFourAnts[3];

   //     -  "[28]" LLN - LLT
   theSum->firstAnt[28]=nadirAnts[0];
   theSum->secondAnt[28]=nextFourAnts[0];

//     -  "[29]" LLN - LT
   theSum->firstAnt[29]=nadirAnts[0];
   theSum->secondAnt[29]=sixAnts[0];
//     -  "[30]" CN - CT
   theSum->firstAnt[30]=nadirAnts[4];
   theSum->secondAnt[30]=sixAnts[1];
//     -  "[31]" CN - LT
   theSum->firstAnt[31]=nadirAnts[4];
   theSum->secondAnt[31]=sixAnts[0];
//     -  "[32]" CN - RT
   theSum->firstAnt[32]=nadirAnts[4];
   theSum->secondAnt[32]=sixAnts[2];
//     -  "[33]" RRN - RT
   theSum->firstAnt[33]=nadirAnts[1];
   theSum->secondAnt[33]=sixAnts[2];
//     -  "[34]" RRN - RRT
   theSum->firstAnt[34]=nadirAnts[1];
   theSum->secondAnt[34]=nextFourAnts[1];
   

   }else{

      //      std::cout << "Here2\t" << nadirAnts[2] << "\t" << nadirAnts[3] << "\n";
//   -  "[19]" LN - RN
   theSum->firstAnt[19]=nadirAnts[2];
   theSum->secondAnt[19]=nadirAnts[3];

   //   -  "[20]" LN - LLB
   theSum->firstAnt[20]=nadirAnts[2];
   theSum->secondAnt[20]=nextFourAnts[2];

  //   -  "[21]" LN - LB
   theSum->firstAnt[21]=nadirAnts[2];
   theSum->secondAnt[21]=sixAnts[3];

  //  -  "[22]" LN - CB
   theSum->firstAnt[22]=nadirAnts[2];
   theSum->secondAnt[22]=sixAnts[4];

//   -  "[23]" LN - RB
   theSum->firstAnt[23]=nadirAnts[2];
   theSum->secondAnt[23]=sixAnts[5];

 //   -  "[24]" RN - RRB
   theSum->firstAnt[24]=nadirAnts[3];
   theSum->secondAnt[24]=nextFourAnts[3];

//   -  "[25]" RN - RB
   theSum->firstAnt[25]=nadirAnts[3];
   theSum->secondAnt[25]=sixAnts[5];

 //   -  "[26]" RN - CB
   theSum->firstAnt[26]=nadirAnts[3];
   theSum->secondAnt[26]=sixAnts[4];

 //   -  "[27]" RN - LB
   theSum->firstAnt[27]=nadirAnts[3];
   theSum->secondAnt[27]=sixAnts[3];


//     -  "[28]" LN - LLT
   theSum->firstAnt[28]=nadirAnts[2];
   theSum->secondAnt[28]=nextFourAnts[0];

//     -  "[29]" LN - LT
   theSum->firstAnt[29]=nadirAnts[2];
   theSum->secondAnt[29]=sixAnts[0];

//     -  "[30]" LN - CT
   theSum->firstAnt[30]=nadirAnts[2];
   theSum->secondAnt[30]=sixAnts[1];

//     -  "[31]" LN - RT
   theSum->firstAnt[31]=nadirAnts[2];
   theSum->secondAnt[31]=sixAnts[2];

//     -  "[32]" RN - RRT
   theSum->firstAnt[32]=nadirAnts[3];
   theSum->secondAnt[32]=nextFourAnts[1];
//     -  "[33]" RN - RT
   theSum->firstAnt[33]=nadirAnts[3];
   theSum->secondAnt[33]=sixAnts[2];
//     -  "[34]" RN - CT
   theSum->firstAnt[34]=nadirAnts[3];
   theSum->secondAnt[34]=sixAnts[1];


   }   
   


   //Now can make correlations and find max, rms, etc.
   for(int corInd=0;corInd<35;corInd++) {
      TGraph *grCor;
      //      std::cout << corInd << "\t" << theSum->firstAnt[corInd] << "\t" << theSum->secondAnt[corInd] << "\n";
      Int_t ci1=AnitaGeomTool::getChanIndexFromAntPol(theSum->firstAnt[corInd],pol);
      Int_t ci2=AnitaGeomTool::getChanIndexFromAntPol(theSum->secondAnt[corInd],pol);

      //      std::cout << nadirAnts[0] << "\t" << corInd << "\t"<< ci1 << " " << ci2 << "  " << theSum->firstAnt[corInd] << "\t" << theSum->secondAnt[corInd] <<std::endl;

      //      if (ci1*ci2<0) continue; // Linda added this condition

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
      int maxIndex=TMath::LocMax(numPoints,theValues);
//       double maxVal=theValues[maxIndex];;
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




CorrelationSummaryAnita3 *PrettyAnitaEvent::createCorrelationSummaryAnita3(Int_t centreAnt,AnitaPol::AnitaPol_t pol, Double_t deltaT)
{
  //Gets the 11 correlations and then takes the max, rms and neighbouring maxima
  if(centreAnt<0)
    centreAnt=getMaxAntenna(pol);

  // Anita 3 has 15 antennas in 5 phi sectors
  // LLT LT CT RT RRT
  // LLM LM CM RM RRM
  // LLB LB CB RB RRB
  // Resulting in 48 correlations

  int nineAnts[9];
  fillNineAntArrays(centreAnt,nineAnts);
  int nextSixAnts[6];
  fillNextSixAntArrays(centreAnt,nextSixAnts);

   CorrelationSummaryAnita3 *theSum = new CorrelationSummaryAnita3(eventNumber, centreAnt, nineAnts,deltaT);
   for(int i=0;i<6;i++)
      theSum->nextSixAnts[i]=nextSixAnts[i];


   //Now need to make correlation index pairs
   //Top-Bottom first
   theSum->firstAnt[0]=nineAnts[0];
   theSum->secondAnt[0]=nineAnts[3];
   theSum->firstAnt[1]=nineAnts[1];
   theSum->secondAnt[1]=nineAnts[4];
   theSum->firstAnt[2]=nineAnts[2];
   theSum->secondAnt[2]=nineAnts[5];
   theSum->firstAnt[3]=nineAnts[3];
   theSum->secondAnt[3]=nineAnts[6];
   theSum->firstAnt[4]=nineAnts[4];
   theSum->secondAnt[4]=nineAnts[7];
   theSum->firstAnt[5]=nineAnts[5];
   theSum->secondAnt[5]=nineAnts[8];
   //Now Left-Right
   theSum->firstAnt[6]=nineAnts[0];
   theSum->secondAnt[6]=nineAnts[1];
   theSum->firstAnt[7]=nineAnts[1];
   theSum->secondAnt[7]=nineAnts[2];
   theSum->firstAnt[8]=nineAnts[3];
   theSum->secondAnt[8]=nineAnts[4];
   theSum->firstAnt[9]=nineAnts[4];
   theSum->secondAnt[9]=nineAnts[5];
   theSum->firstAnt[10]=nineAnts[6];
   theSum->secondAnt[10]=nineAnts[7];
   theSum->firstAnt[11]=nineAnts[7];
   theSum->secondAnt[11]=nineAnts[8];
   //Now Diagonal
   theSum->firstAnt[12]=nineAnts[0];
   theSum->secondAnt[12]=nineAnts[4];
   theSum->firstAnt[13]=nineAnts[2];
   theSum->secondAnt[13]=nineAnts[4];
   theSum->firstAnt[14]=nineAnts[6];
   theSum->secondAnt[14]=nineAnts[4];
   theSum->firstAnt[15]=nineAnts[8];
   theSum->secondAnt[15]=nineAnts[4];
   theSum->firstAnt[16]=nineAnts[1];
   theSum->secondAnt[16]=nineAnts[3];
   theSum->firstAnt[17]=nineAnts[1];
   theSum->secondAnt[17]=nineAnts[5];
   theSum->firstAnt[18]=nineAnts[3];
   theSum->secondAnt[18]=nineAnts[7];
   theSum->firstAnt[19]=nineAnts[7];
   theSum->secondAnt[19]=nineAnts[5];
   //Now Leftmost - centre or centre - rightmost 
   theSum->firstAnt[20]=nextSixAnts[0];
   theSum->secondAnt[20]=nineAnts[1];
   theSum->firstAnt[21]=nineAnts[1];
   theSum->secondAnt[21]=nextSixAnts[1];
   theSum->firstAnt[22]=nextSixAnts[2];
   theSum->secondAnt[22]=nineAnts[4];
   theSum->firstAnt[23]=nineAnts[4];
   theSum->secondAnt[23]=nextSixAnts[3];
   theSum->firstAnt[24]=nextSixAnts[4];
   theSum->secondAnt[24]=nineAnts[7];
   theSum->firstAnt[25]=nineAnts[7];
   theSum->secondAnt[25]=nextSixAnts[5];
   //Now Leftmost - left or right-rightmost
   theSum->firstAnt[26]=nextSixAnts[0];
   theSum->secondAnt[26]=nineAnts[0];
   theSum->firstAnt[27]=nineAnts[2];
   theSum->secondAnt[27]=nextSixAnts[1];
   theSum->firstAnt[28]=nextSixAnts[2];
   theSum->secondAnt[28]=nineAnts[3];
   theSum->firstAnt[29]=nineAnts[5];
   theSum->secondAnt[29]=nextSixAnts[3];
   theSum->firstAnt[30]=nextSixAnts[4];
   theSum->secondAnt[30]=nineAnts[6];
   theSum->firstAnt[31]=nineAnts[8];
   theSum->secondAnt[31]=nextSixAnts[5];
   // Now diagonals to LM and RM
   theSum->firstAnt[32]=nextSixAnts[0];
   theSum->secondAnt[32]=nineAnts[3];
   theSum->firstAnt[33]=nextSixAnts[4];
   theSum->secondAnt[33]=nineAnts[3];
   theSum->firstAnt[34]=nextSixAnts[1];
   theSum->secondAnt[34]=nineAnts[5];
   theSum->firstAnt[35]=nextSixAnts[5];
   theSum->secondAnt[35]=nineAnts[5];
   // Now bottom to top
   theSum->firstAnt[36]=nextSixAnts[0];
   theSum->secondAnt[36]=nextSixAnts[4];
   theSum->firstAnt[37]=nineAnts[0];
   theSum->secondAnt[37]=nineAnts[6];
   theSum->firstAnt[38]=nineAnts[1];
   theSum->secondAnt[38]=nineAnts[7];
   theSum->firstAnt[39]=nineAnts[2];
   theSum->secondAnt[39]=nineAnts[8];
   theSum->firstAnt[40]=nextSixAnts[1];
   theSum->secondAnt[40]=nextSixAnts[5];
   // Now top to bottom displaced right
   theSum->firstAnt[41]=nextSixAnts[0];
   theSum->secondAnt[41]=nineAnts[6];
   theSum->firstAnt[42]=nineAnts[0];
   theSum->secondAnt[42]=nineAnts[7];
   theSum->firstAnt[43]=nineAnts[1];
   theSum->secondAnt[43]=nineAnts[8];
   theSum->firstAnt[44]=nineAnts[2];
   theSum->secondAnt[44]=nextSixAnts[5];
   // Now top to bottom displaced left
   theSum->firstAnt[45]=nineAnts[0];
   theSum->secondAnt[45]=nextSixAnts[4];
   theSum->firstAnt[46]=nineAnts[1];
   theSum->secondAnt[46]=nineAnts[6];
   theSum->firstAnt[47]=nineAnts[2];
   theSum->secondAnt[47]=nineAnts[7];
   theSum->firstAnt[48]=nextSixAnts[1];
   theSum->secondAnt[48]=nineAnts[8];
   // horizontal across center
   theSum->firstAnt[49]  = nineAnts[0];
   theSum->secondAnt[49] = nineAnts[2];
   theSum->firstAnt[50]  = nineAnts[3];
   theSum->secondAnt[50] = nineAnts[5];
   theSum->firstAnt[51]  = nineAnts[6];
   theSum->secondAnt[51] = nineAnts[8];

   return theSum;
}


CorrelationSummaryAnita3 *PrettyAnitaEvent::getCorrelationSummaryAnita3(Int_t centreAnt,AnitaPol::AnitaPol_t pol, Double_t deltaT){

  CorrelationSummaryAnita3 *theSum = createCorrelationSummaryAnita3(centreAnt, pol, deltaT);

   //Now can make correlations and find max, rms, etc.
   for(int corInd=0;corInd<NUM_CORRELATIONS_ANITA3;corInd++) {
      TGraph *grCor;
      //      std::cout << corInd << "\t" << theSum->firstAnt[corInd] << "\t" << theSum->secondAnt[corInd] << "\n";
      Int_t ci1=AnitaGeomTool::getChanIndexFromAntPol(theSum->firstAnt[corInd],pol);
      Int_t ci2=AnitaGeomTool::getChanIndexFromAntPol(theSum->secondAnt[corInd],pol);

      //      std::cout << corInd << "\t"<< ci1 << " " << ci2 << "  " << theSum->firstAnt[corInd] << "\t" << theSum->secondAnt[corInd] <<std::endl;

//       if (ci1*ci2<0) continue; // Linda added this condition

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
      int maxIndex=TMath::LocMax(numPoints,theValues);
//       double maxVal=theValues[maxIndex];

//           Double_t maxVal,rmsVal;
//           Int_t maxIndex;
//           FFTtools::getPeakRmsRectified(grCor,maxVal,rmsVal,&maxIndex);

//           FFTtools::getPeakRmsSqVal(grCor,maxVal,rmsVal,&maxIndex);
//           for(int i=0;i<grCor->GetN();i++) {
// 	    //     	std::cout << i << "\t" << theTimes[i] << "\t" << theValues[i] << "\n";
//      	 if(theValues[i]>maxVal) {
//      	    maxVal=theValues[i];
//      	    maxIndex=i;
//      	 }
//           }
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
   gMinuit->SetFCN(CorSumFCNanita3);
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
//    std::cout << fmin << "\t" << fedm << "\t" << npari << "\t" << nparx 
//    	     << "\t" << istat << std::endl;
   theSum->setFitResults(phiWave,thetaWave,phiWaveErr,thetaWaveErr,fmin);



   return theSum;
}


/*

#define NUM_BINS_PHI 10
#define NUM_BINS_THETA 10
#define NUM_ANTS 2

TCanvas *PrettyAnitaEvent::getTwoAntMap(int ant1,int ant2){

  gStyle->SetPalette(1);

  //Double_t inch=0.0254;//inch in m
  //antenna centre locations from anita note 345.  The coordinates noted are radius to centre of antenna (centre of base) - NEED TO ADJUST FOR PHASE CENTRE ACTUALLY BEING 17.5 CM INTO ANTENNA - and azimuth (az) angle. uses inches and degrees... only upper ring at the moment
  //Double_t antennaCentres[16][2]={{56.831,-45.100},{56.77,-0.363},{56.464,44.59},{56.164,89.801},{56.012,135.095},{56.387,-179.576},{56.59,-134.46},{56.884,-89.728},{50.249,-67.471},{50.218,-22.571},{50.126,22.239},{49.823,67.286},{49.61,112.537},{49.667,157.928},{49.973,-157.101},{50.219,-112.091}};
  //for(int i=0;i<16;i++){
  //  antennaCentres[i][1]+=0.363;//so that everything is calibrated wrt antenna 2
  //  antennaCentres[i][0]/=inch;
  //  antennaCentres[i][0]-=17.5;//converts to phase centre from antenna centre
  //  antennaCentres[i][1]=antennaCentres[i][1]*PI/180; 
  //}

  Double_t PI=3.14159;

  //setup the phi and theta arrays, phi from -pi to pi, theta from -pi/2 to pi/2
  Double_t phiArray[NUM_BINS_PHI];
  Double_t thetaArray[NUM_BINS_THETA];
  for(UInt_t i=0;i<NUM_BINS_PHI;i++){
    phiArray[i] = (i+1/2) * 2*PI/NUM_BINS_PHI - PI;
  }
  for(UInt_t i=0;i<NUM_BINS_PHI;i++){
    thetaArray[i] =(i+1/2) * PI/NUM_BINS_PHI - PI/2;
  }

  //time relative to phi sector 1 (ant 2 and 20)
  Double_t timeRel[NUM_BINS_PHI][NUM_BINS_THETA][NUM_ANTS];
  Double_t gainRel[NUM_BINS_PHI][NUM_BINS_THETA][NUM_ANTS];

  for(int phi=0;phi<NUM_BINS_PHI;phi++){
    for(int theta=0;theta<NUM_BINS_THETA;theta++){

      timeRel[phi][theta] = UsefulAdu5Pat::getDeltaTExpected(ant1,ant2,phiArray[phi],thetaArray[theta],10.);

    }
  }

  Double_t deltaT=1/(2.6*8);

  TGraph *gr[NUM_ANTS]={0};

  gr[0] = getInterpolatedGraph(AnitaGeomTool::getChanFromAntPol(ant1,AnitaPol::kVertical),deltaT);
  gr[1] = getInterpolatedGraph(AnitaGeomTool::getChanFromAntPol(ant2,AnitaPol::kVertical),deltaT);

  Int_t numBins=gr[0]->GetN();

  Double_t maxPeak[NUM_BINS_PHI][NUM_BINS_THETA];
  Double_t testPeak;
  Double_t testX[NUM_ANTS]={0};
  Double_t testY[NUM_ANTS]={0};

  for(int phi=0;phi<NUM_BINS_PHI;phi++){
    for(int theta=0;theta<NUM_BINS_THETA;theta++){

      maxPeak[phi][theta]=0.;
      testPeak=0.;

      Int_t binShift=static_cast<Int_t>((timeRel[phi][theta][ant2]-timeRel[phi][theta][ant1])/deltaT);

      for(int bin=0;bin<numBins;bin++){

	if((bin+binShift)<0) continue;
	if((bin+binShift)>numBins) continue;

	gr[0]->GetPoint(bin,testX[0],testY[0]);
	gr[1]->GetPoint(bin+binShift,testX[1],testY[1]);

	testPeak=testY[0]+testY[1];
	if(testPeak>maxPeak[phi][theta]) maxPeak[phi][theta]=testPeak;


      }
     
    }
  }

  TCanvas *canvas = (TCanvas*) gROOT->FindObject("twoAntCorr");
  if(!canvas) {
    canvas = new TCanvas("twoAntCorr","twoAntCorr",600,500);
  }

  TH2F *corrPlot = new TH2F("Two Ant Correlation","Two Ant Correlation",NUM_BINS_PHI,0,2*PI,NUM_BINS_THETA,-PI/2,PI/2);

  for(int phi=0;phi<NUM_BINS_PHI;phi++){
    for(int theta=0;theta<NUM_BINS_THETA;theta++){
      corrPlot->Fill(phiArray[phi],thetaArray[theta],maxPeak[phi][theta]);
    }
  }

  corrPlot->Draw("colz");

  return canvas;

}
*/
