//////////////////////////////////////////////////////////////////////////////
/////  PrettyAnitaEvent.h        Useful ANITA event class                      /////
/////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for plotting stuff like waveforms and correlations/////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#ifndef PRETTYANITAEVENT_H
#define PRETTYANITAEVENT_H
#include "UsefulAnitaEvent.h"
#include "CorrelationSummary.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <string>

class PrettyAnitaEvent: public UsefulAnitaEvent
{


 public:
  PrettyAnitaEvent(RawAnitaEvent *eventPtr,WaveCalType::WaveCalType_t calType, PrettyAnitaHk *theHk);


  //Putative Analysis methods
  int getMaxAntenna(AnitaPol::AnitaPol_t pol);
  CorrelationSummary *getCorrelationSummary(Int_t centreAnt,AnitaPol::AnitaPol_t pol,Double_t deltaT=0);


  //Canvas panel getters
  TCanvas *getSixWaveformCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  TCanvas *getTenWaveformCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);

  TCanvas *getSixFFTPowerCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  TCanvas *getSixPowerEnvelopeCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  TCanvas *getSixInterpolatedCanvas(int ant, AnitaPol::AnitaPol_t pol, Double_t deltaT=(1./(2.6*8)), TCanvas *can=0);
  TCanvas *getSixCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  TCanvas *getTenCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  TCanvas *getElevenCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  TCanvas *getElevenInterpolationCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol,Double_t deltaT=(1./(2.6*8)), TCanvas *can=0);
  TCanvas *getSixInterpolatedCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol,Double_t deltaT=(1./(2.6*8)), TCanvas *can=0);


  //Graph getters
  TGraph *getSimplePowerEnvelopeGraph(int chanIndex);
  TGraph *getInterpolatedGraph(int chanIndex, double deltaT);
  TGraph *getFFTMagnitude(TGraph *grIn); 
  TGraph *getFFTMagnitude(int chanIndex); 
  TGraph *getCorrelation(int chanIndex1, int chanIndex2); 
  TGraph *getCorrelation(TGraph *gr1, TGraph *gr2); 
  TGraph *getCorrelationInterpolated(int chanIndex1, int chanIndex2, Double_t deltaT=(1./(2.6*8)) ); 
  




  ClassDef(PrettyAnitaEvent,1);

 private:
  void setStyleSixCanvas();
  void fillSixAntArrays(int ant, int topAnts[3], int bottomAnts[3]);
  void fillNextFourAntArrays(int ant, int nextFourAnts[4]);
  int getPrettyColour(int index);

  Double_t fDeltaT;
  Double_t fWaveOffset;


};


#endif //PRETTYANITAEVENT_H
