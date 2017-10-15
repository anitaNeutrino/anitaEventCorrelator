/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
 Class to greatly simplify the drawing of pretty TGraph of Antarctica.
 It is intended as a partial replacement to AntarcticaMapPlotter in anitaAnalysisTools.
 Uses the new and improved RampdemReader to do the hard coordinate conversion.
***********************************************************************************************************/

#ifndef TGRAPHANTARCTICA_H
#define TGRAPHANTARCTICA_H

#include "TGraph.h"
#include "RampdemReader.h"
#include "TGraphEditor.h"
#include "AntarcticaBackground.h"
#include "TChain.h"
#include "TCut.h"
#include "BaseList.h"


/**
 * @class TGraphAntarctica
 *
 * Class to greatly simplify the task of drawing TGraphs on a background of Antarctica
 * This class expects you to provide lon, lat as the X, Y inputs. Get it the right way round!
 * It then talks behind the scenes to RampdemReader and converts those numbers to Easting/Northing
 * The Draw option has been redefined to plot everything on a background of Antarctica.
 *
 * Getting something approximating standard TGraph drawing behaviour was a massive pain...
 * There's probably a few bugs in here still, but I'm happy enough with the current functionality.
 *
 */
class TGraphAntarctica : public TGraph {


public:

  enum {
    defaultGpsTreeStride = 10000
  };

  // static members
  static TGraphAntarctica* makeGpsGraph(int firstRun, int lastRun, int gpsTreeStride=defaultGpsTreeStride, bool quiet=true);

  // boring constructors
  TGraphAntarctica() : TGraph() {init();}
  explicit TGraphAntarctica(Int_t n) : TGraph(n) {init();}
  TGraphAntarctica(Int_t n, const Int_t *x, const Int_t *y) : TGraph(n, x, y) {init();}
  TGraphAntarctica(Int_t n, const Float_t *x, const Float_t *y) : TGraph(n, x, y) {init();}
  TGraphAntarctica(Int_t n, const Double_t *x, const Double_t *y) : TGraph(n, x, y) {init();}
  explicit TGraphAntarctica(const TGraph &gr) : TGraph(gr) {init();}

  // interesting constructors
  TGraphAntarctica(TChain* chain, TString lonSelector="longitude", TString latSelector="latitude", TCut cut = "");
  TGraphAntarctica(TTree* tree, TString lonSelector="longitude", TString latSelector="latitude", TCut cut = "");

  explicit TGraphAntarctica(const TVectorF &vx, const TVectorF &vy) : TGraph(vx, vy) {init();}

  explicit TGraphAntarctica(const TVectorD &vx, const TVectorD &vy) : TGraph(vx, vy) {init();}
  explicit TGraphAntarctica(const TH1 *h) : TGraph(h){init();}

  explicit TGraphAntarctica(const TF1 *f, Option_t *option="") : TGraph(f, option){init();}

  explicit TGraphAntarctica(const char *filename, const char *format="%lg %lg", Option_t *option="") : TGraph(filename, format, option) {init();}
  TGraphAntarctica(const BaseList::base& b);


  virtual void SetPoint(Int_t i, Double_t lon, Double_t lat);

  Double_t* GetEasting(){return GetX();}
  Double_t* GetNorthing(){return GetY();}
  void Draw(Option_t* opt = "");
  void ExecuteEvent(Int_t event, Int_t x, Int_t y);

  ClassDef(TGraphAntarctica,1)

private:

  void init();
  Bool_t doneConversion; // convert initial array from lat/lon to easting northing
  void convertArrays();

  TString fToolTipExtraTextFormat;

};







#endif
