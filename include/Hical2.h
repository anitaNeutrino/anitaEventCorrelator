/*
getting rid of hical 2.
 */

#ifndef HC2_H
#define HC2_H

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "HCHKTree.h"
#include "AnitaDataset.h"
#include "AnitaVersion.h"
#include "FFTtools.h"
#include "AnitaEventSummary.h"
#include "RawAnitaHeader.h"
#include "Adu5Pat.h"

#include <iostream>

using namespace std;

class Hical2{
public:

  Hical2();
  //  ~Hical2();
  //the main cut function using event number
  double isHical(UInt_t eventNumber, int geomCut=1);
  //using event summary
  static  double isHical(AnitaEventSummary *sum, int geomCut=1);
  //will return dphi for specified payload for map peaks 1, 2, 3
  double dPhi(int aorb, int peak=0);
  //return the absolute angles to HC
  static int whereAreHical(UInt_t eventNumber, double * angleToA, double *angleToB);
  //is hical on
  static bool hc2aOn(UInt_t triggerTime);
  static bool hc2bOn(UInt_t triggerTime);
  
private:

  //called on the first call of isHical, loads the root files into the trees
  int initHical();
  int INIT_HICAL=0;

  //angle things
  static double angleToThing(double lat1, double lon1, double lat2, double lon2);
  static double deg2rad(double deg);
  static double rad2deg(double rad);

  //anitadataset, and trees, and storage stuff to make things faster
  AnitaDataset *adset =0;
  //TTree * hc2ahk_tree = new TTree();
  //TTree * hc2bhk_tree = new TTree();
  
  //HCHKTree *hc2ahk=new HCHKTree();
  //HCHKTree *hc2bhk=new HCHKTree();


  ClassDef(Hical2, 1);
};
#endif
