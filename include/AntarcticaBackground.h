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

// I want a thing in the background I can click on and send commands to -> must be a TH2D
// Easy to pick biggest thing and read in data as changed...
// What's hard is to rebin the internals dynamically


class AntarcticaBackground : public TProfile2D {

public:


  AntarcticaBackground();  ///< Default constructor... don't use this

  // sensible constructor, but still don't use this
  AntarcticaBackground(RampdemReader::dataSet dataSet, Int_t coarseness,
		       const char *name,const char *title,
		       Int_t nx, Double_t xlow, Double_t xup,
		       Int_t ny, Double_t ylow, Double_t yup,
		       Option_t *option = "");

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

  // void makePrettyPalette();
  // void setPadMargins();
  // void setColAxisTitle();

  void update();

public:

  // GUI Stuff... For the context menus


  // Here I implement the context menu functions...
  // there should be one getter/setter boolian (*TOGGLE*) function

  // I'm giving the option of the enabled data sets in RampdemReader.
  // rampdem, bed, icemask_grounded_and_shelves, surface, thickness,

  void SetRampdemDataSet(bool useRampdemDataSet); //*TOGGLE* *GETTER=GetRampdemDataSet
  Bool_t GetRampdemDataSet();

  void SetBedDataSet(bool useBedDataSet); //*TOGGLE* *GETTER=GetBedDataSet
  Bool_t GetBedDataSet();

  void SetIcemaskDataSet(bool useIcemaskDataSet); //*TOGGLE* *GETTER=GetIcemaskDataSet
  Bool_t GetIcemaskDataSet();

  void SetSurfaceDataSet(bool useSurfaceDataSet); //*TOGGLE* *GETTER=GetSurfaceDataSet
  Bool_t GetSurfaceDataSet();

  void SetThicknessDataSet(bool useThicknessDataSet); //*TOGGLE* *GETTER=GetThicknessDataSet
  Bool_t GetThicknessDataSet();

  void SetDrawLonLatGrids(Bool_t drawLonLatGrids); //*TOGGLE* *GETTER=GetDrawLonLatGrids
  Bool_t GetDrawLonLatGrids();

  Int_t GetCoarseness();
  void SetCoarseness(Int_t coarseness); // *MENU* *ARGS={coarseness=>fCoarseness}

  RampdemReader::dataSet GetDataSet();
  void SetDataSet(RampdemReader::dataSet dataSet);

  static AntarcticaBackground* generate(RampdemReader::dataSet dataSet, Int_t coarseness); // use this instead of constructor

  ClassDef(AntarcticaBackground, 1)


};



#endif
