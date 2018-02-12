/*
finding/getting rid of hical 2.
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

#include <cmath>
#include <iostream>

//#define pi TMath::Pi()
using namespace std;

class Hical2{
public:

  Hical2();
  //  ~Hical2();
  
  //cut function using event number, slow but thourough, using pointing, H/V ratio, and template matching
  double isHical(UInt_t eventNumber);
  
  //using event summary. a purely 5 sigma pointing cut on the indicated HPol peak.
  static double isHical(AnitaEventSummary *sum, int peak=0);
  //to be used inside of TTree::Draw, giving it an event number and absolute direction for the peak
  static double isHical(UInt_t eventNumber, double peakPhi);
  //to be used inside of TTree::Draw, giving it an event number, trigger time, and absolute direction for the peak
  static double isHical(UInt_t eventNumber, UInt_t triggerTime, double peakPhi);
    //to be used inside of TTree::Draw, giving it an event number, trigger time, and absolute direction for the peak
  static double isHical(UInt_t eventNumber, UInt_t triggerTime, UInt_t triggerTimeNs, double peakPhi);
  
  /* //will return dphi for specified payload for map peaks 1, 2, 3 */
  /* double dPhi(int aorb, int peak=0); */
  
  //find lat, lon, alt for both payloads
  static int whereAreHical(UInt_t eventNumber, double *latA, double *lonA, double *altA, double *latB, double *lonB, double *altB);
  //same function as above
  static int whereIsHical(UInt_t eventNumber, double *latA, double *lonA, double *altA, double *latB, double *lonB, double *altB);
  //find the absolute angles to HC  
  static int angleToHical(UInt_t eventNumber, double * angleToA, double *angleToB);
  //same function in case you forget the plurality of hical
  static int anglesToHical(UInt_t eventNumber, double * angleToA, double *angleToB);
  //is hical on
  static bool hc2aOn(UInt_t triggerTime);
  static bool hc2bOn(UInt_t triggerTime);
  
private:

  //called on the first call of isHical, loads the root files into the trees
  int initHical();
  int INIT_HICAL=0;

  //utility functions
  static double angleToThing(double lat1, double lon1, double lat2, double lon2);
  static double straightLineDist(double lat1, double lon1, double alt1, double lat2, double lon2, double alt2);
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
