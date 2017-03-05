/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
 TProfile2D to contain the map of antarctica for pretty, interactive plotting.
 Honed to "perfection" (well, a functional state) on a flight from London to Los Angeles.
***********************************************************************************************************/

#ifndef ANTARCTICA_BACKGROUND_H
#define ANTARCTICA_BACKGROUND_H

#include "TProfile2D.h"
#include "RampdemReader.h"
#include "TGraph.h"
#include "TGToolTip.h"

const int defaultCoarseness = 10; // 1 is v slow on my laptop but you might want to do that when you zoom in.

class AntarcticaBackground : public TProfile2D {

public:

  AntarcticaBackground(RampdemReader::dataSet dataSet = RampdemReader::thickness,
		       Int_t coarseness = defaultCoarseness);
  void Draw(Option_t* opt = "colz");


  // Interactive plotting fun
  // Click on the map to choose which RAMPDEM/BEDMAP2 data and set resolution.
  Int_t GetCoarseness();
  void SetCoarseness(Int_t coarseness); // *MENU* *ARGS={coarseness=>fCoarseness}

  void Rampdem(bool useRampdem); //*TOGGLE* *GETTER=GetRampdem
  Bool_t GetRampdem();

  void Bed(bool useBed); //*TOGGLE* *GETTER=GetBed
  Bool_t GetBed();

  void Icemask(bool useIcemask); //*TOGGLE* *GETTER=GetIcemask
  Bool_t GetIcemask();

  void Surface(bool useSurface); //*TOGGLE* *GETTER=GetSurface
  Bool_t GetSurface();

  void Thickness(bool useThickness); //*TOGGLE* *GETTER=GetThickness
  Bool_t GetThickness();

  void Grid(Bool_t grid); //*TOGGLE* *GETTER=GetGrid
  Bool_t GetGrid();

  void SetGridDivisions(Int_t deltaLon, Int_t deltaLat); // *MENU* *ARGS={deltaLat=>fDeltaLon, deltaLon=>fDeltaLat}

  void ToolTip(Bool_t toolTip); // *TOGGLE* *GETTER=GetToolTip
  Bool_t GetToolTip();


  RampdemReader::dataSet GetDataSet();
  void SetDataSet(RampdemReader::dataSet dataSet);
  // void Interactive(Int_t event, Int_t x, Int_t y, TObject* selected);
  void ExecuteEvent(Int_t event, Int_t x, Int_t y); //, TObject* selected);

private:

  Int_t fCoarseness; // map coarseness factor
  RampdemReader::dataSet fDataSet; // rampdem data set
  Bool_t needRemakeHist; // if changed coarseness or data set

  Bool_t fGrid;
  Int_t fDeltaLon;
  Int_t fDeltaLat;
  void updateGrid();
  void deleteGrid();
  std::vector<TGraph*> grGrids; ///< The internally stored lat/lon grids.

  Int_t fGridPoints; ///< Number of points on the lat/lon grid.
  Bool_t needRemakeGrid; // if changed grid settings
  void updateHist();
  void init(RampdemReader::dataSet dataSet, Int_t coarseness); // in case I ever want another constructor

  void setPadMargins(); // prettification
  void prettifyPalette(); // prettification

  Bool_t fAlreadyDrawn; // Set by Draw(), help the updateHist() function to do sensible things
  ClassDef(AntarcticaBackground, 0)

  Bool_t fUseToolTip;
  TGToolTip* fToolTip;
  TString fToolTipUnits;
  void setToolTipUnits();
};



#endif
