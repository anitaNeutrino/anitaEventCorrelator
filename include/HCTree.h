#ifndef HC2_TREE
#define HC2_TREE

#include "TGraph.h"

class HCTree {

 public:

  HCTree();

  bool on;//is hc on?

  UInt_t run;//a4 run
  UInt_t eventNumber;//a4 event number
  UInt_t triggerTime;//a4 time in seconds
  UInt_t triggerTimeNs;//a4 time in nanoseconds

  double a4TimeD;//a4 time formatted to be seconds.subseconds
  
  
  UInt_t dTimeHk;//time difference between a4 timestamp and hical HK timestamp
  double dTimeSci;//time difference between a4 timestamp and hical science (pulse) timestamp
  UInt_t dTimeMip;//time difference between a4 timestamp and the hical mip data
  
  //everything below depends upon the above two variables, as the following variables are not interpolated values. these are filled with the hical values at the smallest d_time_hk or d_time_sci (depending upon what variable is being filled).

  UInt_t hcTime;//hc time in seconds
  UInt_t hcTimeNs;//hc time in nanoseconds
  double hcTimeD;//hc time as a double  

  double dt;//anita time - hical time - time of flight
  
  double lat;
  double lon;
  double alt;
  double az;//absolute orientation of HC antanna w/r/t anita
  double a4_hc;//angle from anita to hical
  double hc_a4;//angle from hical to anita
  
  double dist;//straight line-of-sight distance
  double distLatLon;//ground distance between payloads
  double thetaD;//calculated direct angle to HC
  double thetaR;//calculated angle to reflected event were there to be one
  double thetaI;

  //zeroth peak
  double dThetaD;//the difference between the anita peak theta and hc's theta D
  double absDThetaD;
  double dThetaR;//the difference between the anita peak theta and HC's theta R
  double absDThetaR;
  double dPhi;//the difference between the anita peak bearing and hc's bpsbearing
  double absDPhi;
  //first peak
  double dThetaD1;//the difference between the anita peak theta and hc's theta D 
  double absDThetaD1;
  double dThetaR1;//the difference between the anita peak theta and HC's theta R 
  double absDThetaR1;
  double dPhi1;//the difference between the anita peak bearing and hc's bpsbearing
  double absDPhi1;
  //second peak
  double dThetaD2;//the difference between the anita peak theta and hc's theta D
  double absDThetaD2;
  double dThetaR2;//the difference between the anita peak theta and HC's theta R
  double absDThetaR2;
  double dPhi2;//the difference between the anita peak bearing and hc's bpsbearing
  double absDPhi2;

  double tofDir;
  UInt_t tofDirNs;
  double tofRefl;
  UInt_t tofReflNs;
  double tDiff;
  double tDiffNs;

  double a4TDiff1;//these are the dt of nanosecond timing between this event in a4 and the previous 1...5
  double a4TDiff2;
  double a4TDiff3;
  double a4TDiff4;
  double a4TDiff5;

  UInt_t a4TDiff1Ns;
  UInt_t a4TDiff2Ns;
  UInt_t a4TDiff3Ns;
  UInt_t a4TDiff4Ns;
  UInt_t a4TDiff5Ns;
  
  double volts;//battery volts
  double pressure;//pv pressure, only valid when on.
};

inline HCTree::HCTree(){}

class SimpleHCTree{
public:
  SimpleHCTree();

  UInt_t run=0;
  UInt_t eventNumber=0;
  int on=0;
  
  UInt_t dTimeHk=0;//time difference between a4 timestamp and hical HK timestamp
  double dTimeSci=0;//time difference between a4 timestamp and hical science (pulse) timestamp
  UInt_t dTimeMip=0;//time difference between a4 timestamp and the hical mip data

  double lat=0;
  double lon=0;
  double alt=0;
  
  double a4_hc=0;//angle from anita to hical
  
  double dt=0;
  double dPhi=0;
  double absDPhi=0;
  double dPhi1=0;
  double absDPhi1=0;
  double dPhi2=0;
  double absDPhi2=0;

  double hilbertHVRatio=0;
  //  double peakHilbertHPol;
  //double peakHilbertVPol;
  
};
inline SimpleHCTree::SimpleHCTree(){}

class HCPairTree {

 public:

  HCPairTree();

  UInt_t run;
  UInt_t devtno;
  UInt_t revtno;
  double dtheta;
  double peakdh;
  double peakdv;
  double peakrh;
  double peakrv;
  double rtheta;
  double peakHilbertDH;
  double peakHilbertRH;
  double peakHilbertDV;
  double peakHilbertRV;
  double dphi;
  double rphi;
  double dtime;
  TGraph *d=0;
  TGraph *r=0;
  TGraph *dv=0;
  TGraph *rv=0;
  TGraph *dc=0;
  TGraph *rc=0;
  double az;
  double dtimesci;
  double tdiff;
};

inline HCPairTree::HCPairTree(){}

class CorrTree{
public:

  CorrTree();

  double dhh;
  double dvv;
  double dhv;
  double dvh;
  double rhh;
  double rvv;
  double rhv;
  double rvh;
  double weight;
  double peakh;
  double peakv;
};

inline CorrTree::CorrTree(){}



class HCHKTree{
 public:
  HCHKTree();
  double lat;
  double lon;
  double alt;
  double time;
  //ClassDef(HCHKTree, 0);
};
//ClassImp(HCHKTree)
inline HCHKTree::HCHKTree(){}

#endif
