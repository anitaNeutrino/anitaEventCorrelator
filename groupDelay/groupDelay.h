//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 14 11:16:03 2009 by ROOT version 5.20/00
// from TTree deltaTTree/Tree of Delta T's
// found on file: /Users/rjn/anita/data/deltaTTrees/justTaylorKurtNumbers/deltaTFile13.root
//////////////////////////////////////////////////////////

#ifndef groupDelay_h
#define groupDelay_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "AnitaGeomTool.h"

class groupDelay {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Long64_t        entry;
   Int_t           firstAnt;
   Int_t           secondAnt;
   Int_t           maxAnt;
   Int_t           labChip;
   Double_t        deltaT;
   Double_t        deltaTExpected;
   Double_t        corPeak;
   Double_t        corRMS;
   Double_t        phiMaxAnt;
   Double_t        phiWave;
   Double_t        thetaWave;
   UInt_t          eventNumber;
   UInt_t          triggerTime;
   UInt_t          triggerTimeNs;
   Int_t           corInd;
   Double_t        balloonLat;
   Double_t        balloonLon;
   Double_t        balloonAlt;
   Double_t        heading;
   Double_t        pitch;
   Double_t        roll;
   Double_t        meanPhiAntPair;
   Double_t        deltaPhiAntPair;
   UInt_t          expTaylorTime;

   // List of branches
   TBranch        *b_entry;   //!
   TBranch        *b_firstAnt;   //!
   TBranch        *b_secondAnt;   //!
   TBranch        *b_maxAnt;   //!
   TBranch        *b_labChip;   //!
   TBranch        *b_deltaT;   //!
   TBranch        *b_deltaTExpected;   //!
   TBranch        *b_corPeak;   //!
   TBranch        *b_corRMS;   //!
   TBranch        *b_phiMaxAnt;   //!
   TBranch        *b_phiWave;   //!
   TBranch        *b_thetaWave;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_triggerTime;   //!
   TBranch        *b_triggerTimeNs;   //!
   TBranch        *b_corInd;   //!
   TBranch        *b_balloonLat;   //!
   TBranch        *b_balloonLon;   //!
   TBranch        *b_balloonAlt;   //!
   TBranch        *b_heading;   //!
   TBranch        *b_pitch;   //!
   TBranch        *b_roll;   //!
   TBranch        *b_meanPhiAntPair;   //!
   TBranch        *b_deltaPhiAntPair;   //!
   TBranch        *b_expTaylorTime;   //!

   groupDelay(TTree *tree=0);
   virtual ~groupDelay();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef groupDelay_cxx
groupDelay::groupDelay(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/rjn/anita/data/deltaTTrees/justTaylorKurtNumbers/deltaTFile13.root");
      if (!f) {
         f = new TFile("/Users/rjn/anita/data/deltaTTrees/justTaylorKurtNumbers/deltaTFile13.root");
      }
      tree = (TTree*)gDirectory->Get("deltaTTree");

   }
   Init(tree);
}

groupDelay::~groupDelay()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t groupDelay::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t groupDelay::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void groupDelay::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("entry", &entry, &b_entry);
   fChain->SetBranchAddress("firstAnt", &firstAnt, &b_firstAnt);
   fChain->SetBranchAddress("secondAnt", &secondAnt, &b_secondAnt);
   fChain->SetBranchAddress("maxAnt", &maxAnt, &b_maxAnt);
   fChain->SetBranchAddress("labChip", &labChip, &b_labChip);
   fChain->SetBranchAddress("deltaT", &deltaT, &b_deltaT);
   fChain->SetBranchAddress("deltaTExpected", &deltaTExpected, &b_deltaTExpected);
   fChain->SetBranchAddress("corPeak", &corPeak, &b_corPeak);
   fChain->SetBranchAddress("corRMS", &corRMS, &b_corRMS);
   fChain->SetBranchAddress("phiMaxAnt", &phiMaxAnt, &b_phiMaxAnt);
   fChain->SetBranchAddress("phiWave", &phiWave, &b_phiWave);
   fChain->SetBranchAddress("thetaWave", &thetaWave, &b_thetaWave);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("triggerTime", &triggerTime, &b_triggerTime);
   fChain->SetBranchAddress("triggerTimeNs", &triggerTimeNs, &b_triggerTimeNs);
   fChain->SetBranchAddress("corInd", &corInd, &b_corInd);
   fChain->SetBranchAddress("balloonLat", &balloonLat, &b_balloonLat);
   fChain->SetBranchAddress("balloonLon", &balloonLon, &b_balloonLon);
   fChain->SetBranchAddress("balloonAlt", &balloonAlt, &b_balloonAlt);
   fChain->SetBranchAddress("heading", &heading, &b_heading);
   fChain->SetBranchAddress("pitch", &pitch, &b_pitch);
   fChain->SetBranchAddress("roll", &roll, &b_roll);

   fChain->SetBranchAddress("meanPhiAntPair", &meanPhiAntPair, &b_meanPhiAntPair);
   fChain->SetBranchAddress("deltaPhiAntPair", &deltaPhiAntPair, &b_deltaPhiAntPair);
   fChain->SetBranchAddress("expTaylorTime", &expTaylorTime, &b_expTaylorTime);
   Notify();
}

Bool_t groupDelay::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void groupDelay::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t groupDelay::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  if(corRMS<40) return -1;
  if((corPeak/corRMS)<10) return -1;
  //  if(firstAnt>32) return -1;

  int middleAnt=firstAnt; 
  int leftAnt,rightAnt;
  fGeomTool->getThetaPartners(middleAnt,leftAnt,rightAnt); 
  if(secondAnt!=rightAnt) return -1;
  
  return 1;
}
#endif // #ifdef groupDelay_cxx
