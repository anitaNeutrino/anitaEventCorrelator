#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "UsefulAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "UsefulAdu5Pat.h"
#include "CorrelationSummary.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>

void plotFakePrettyThings(int entry, int ant);

  
void plotFakePrettyThings(int entry, int ant) {

  char eventName[FILENAME_MAX];
  sprintf(eventName,"/home/rjn/anita/nutsAndBolts/naughty/peterFakeEvents.root");
  RawAnitaEvent *event = 0;


  TChain *eventChain = new TChain("eventTree");
  eventChain->Add(eventName);
  eventChain->SetBranchAddress("event",&event);


  //Friends only seem to work with TTree::Draw and similar commands
  //if you are manually calling GetEntry (i.e in a loop) you must call
  //the GetEntry for each tree separately.
  //  eventChain->AddFriend(headTree);
  //  eventChain->AddFriend(prettyHkTree);
    
  //Stupidly most do this to be perfectly safe  
  eventChain->GetEntry(entry);
    
  PrettyAnitaHk *hk=0;
  PrettyAnitaEvent realEvent(event,WaveCalType::kVTFullJWPlusClock,hk);
  
  //  for(int samp=0;samp<260;samp++) {
  //    cout << samp << "\t" << realEvent.fTimes[0][samp] << "\t" << realEvent.fVolts[0][samp] << endl;
  //  }

//  TCanvas *canWave = realEvent.getSixWaveformCanvas(ant, AnitaPol::kVertical);
  TCanvas *canWave = realEvent.getTenWaveformCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canPower = realEvent.getSixPowerEnvelopeCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canInter = realEvent.getSixInterpolatedCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canFFT = realEvent.getSixFFTPowerCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canCor = realEvent.getSixCorrelationCanvas(ant, AnitaPol::kVertical);
  TCanvas *canElevenCor = realEvent.getElevenCorrelationCanvas(ant, AnitaPol::kVertical);
  //  TCanvas *canElevenIntCor = realEvent.getElevenInterpolationCorrelationCanvas(ant, AnitaPol::kVertical);
 // TCanvas *canIntCor = realEvent.getSixInterpolatedCorrelationCanvas(ant, AnitaPol::kVertical);

}

