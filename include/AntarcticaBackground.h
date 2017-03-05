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


const int defaultCoarseness = 10; // 1 is v slow on my laptop but you might want to do that when you zoom in.

class AntarcticaBackground : public TProfile2D {

public:

  AntarcticaBackground(RampdemReader::dataSet dataSet = RampdemReader::rampdem,
		       Int_t coarseness = defaultCoarseness);
  void Draw(Option_t* opt = "colz");


  // Interactive plotting fun
  // Click on the map to choose which RAMPDEM/BEDMAP2 data set to plot.

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

  void DrawLonLatGrids(Bool_t drawLonLatGrids); //*TOGGLE* *GETTER=GetDrawLonLatGrids
  Bool_t GetDrawLonLatGrids();

  Int_t GetCoarseness();
  void SetCoarseness(Int_t coarseness); // *MENU* *ARGS={coarseness=>fCoarseness}

  RampdemReader::dataSet GetDataSet();
  void SetDataSet(RampdemReader::dataSet dataSet);

private:

  Int_t fCoarseness; // map coarseness factor
  Int_t lastCoarseness; // remember last drawn coarseness factor

  RampdemReader::dataSet fDataSet; // rampdem data set
  RampdemReader::dataSet lastDataSet; // remember last drawn data set

  // void deleteLonLatGrids();
  // void makeLonLatGrids();
  Bool_t fDrawLonLatGrids;
  // std::vector<TGraph*> grLonLatGrids; ///< The internally stored lat/lon grids.
  // Int_t fLonLatGridPoints; ///< Number of points on the lat/lon grid.
  // Int_t lastLonLatGridPoints; ///< Remember last drawn points for lat/lon grid.
  // void setColAxisTitle();

  void update();
  void init(RampdemReader::dataSet dataSet, Int_t coarseness); // in case I ever want another constructor

  void setPadMargins(); // prettification
  void prettifyPalette(); // prettification


  ClassDef(AntarcticaBackground, 1)


};



#endif
