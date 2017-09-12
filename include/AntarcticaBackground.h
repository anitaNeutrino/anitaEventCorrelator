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
#include "RVersion.h"
#include "TColor.h"

class TExec;
class TGToolTip;
class TGraphAntarctica;
class TPad;

namespace AntarcticaBackgroundDefaults {
const int defaultCoarseness = 10;  /// 1 is v slow on my laptop maybe set lower to zoom in
const double zAxisTextSize = 0.02; /// Standard size for the z-axis text
const double zAxisWidth = 0.03; /// Standard size for z-xaxis
const double zAxisHeight = 0.4; /// Standard size for z-xaxis
const double zAxisRightMargin = 0.02; /// Standard position for z-xaxis
const double zAxisTopBottomMargin = 0.02; /// Standard position for z-xaxis
}


class AntarcticaBackground : public TProfile2D {


 public:

  AntarcticaBackground(RampdemReader::dataSet dataSet = RampdemReader::thickness,
                       Int_t coarseness = AntarcticaBackgroundDefaults::defaultCoarseness);
  virtual ~AntarcticaBackground();

  void Draw(Option_t* opt = "colz");

  // Interactive plotting fun
  // Enable for extra info on mouse over...
  void UnZoom(){fXaxis.UnZoom();fYaxis.UnZoom();} // *MENU

  void SetGrayScale(bool greyScale); //*TOGGLE* *GETTER=GetGrayScale
  Bool_t GetGrayScale() const;

  void SetToolTip(Bool_t toolTip); // *TOGGLE* *GETTER=GetToolTip
  Bool_t GetToolTip() const;

  void SetShowBases(Bool_t showBases); // *TOGGLE* *GETTER=GetShowBases
  Bool_t GetShowBases() const;

  // Click on the map to choose which RAMPDEM/BEDMAP2 data and set resolution.
  Int_t GetCoarseness() const;
  void SetCoarseness(Int_t coarseness); // *MENU* *ARGS={coarseness=>fCoarseness}

  void SetRampdem(bool useRampdem); //*TOGGLE* *GETTER=GetRampdem
  Bool_t GetRampdem() const;

  void SetBed(bool useBed); //*TOGGLE* *GETTER=GetBed
  Bool_t GetBed() const;

  void SetIcemask(bool useIcemask); //*TOGGLE* *GETTER=GetIcemask
  Bool_t GetIcemask() const;

  void SetSurface(bool useSurface); //*TOGGLE* *GETTER=GetSurface
  Bool_t GetSurface() const;

  void SetThickness(bool useThickness); //*TOGGLE* *GETTER=GetThickness
  Bool_t GetThickness() const;

  void SetGrid(Bool_t grid); //*TOGGLE* *GETTER=GetGrid
  Bool_t GetGrid() const;

  void SetGridDivisions(Int_t deltaLon, Int_t deltaLat); // *MENU* *ARGS={deltaLat=>fDeltaLon, deltaLon=>fDeltaLat}
  void PrettifyColorAxis(); // *MENU

  const char* getToolTipUnits(){return fToolTipUnits.Data();}

  RampdemReader::dataSet GetDataSet() const;
  void SetDataSet(RampdemReader::dataSet dataSet);
  void ExecuteEvent(Int_t event, Int_t x, Int_t y);
  void updateToolTip(Int_t event, Int_t x, Int_t y, const char* extraInfo = NULL);

  // needs to be public and accessible for other classes that want to find one of these on their canvases
  static const char* getDefaultName(){return "fAntarctica";}

  void setPalette();
  void unsetPalette();

  // at some point, supporting ROOT versions < 6 is gonna be impossible...
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    std::map<RampdemReader::dataSet, EColorPalette> palettes; //! Does not persist in ROOT!
    float opacity; //! Does not persist in ROOT!
#endif


 private:


    Int_t fCoarseness; // map coarseness factor
    RampdemReader::dataSet fDataSet; // rampdem data set
    Bool_t needRemakeHist; // if changed coarseness or data set
    Bool_t fGrayScale;


    Bool_t fGrid;
    Int_t fDeltaLon;
    Int_t fDeltaLat;
    void updateGrid();
    void deleteGrid();
    std::vector<TGraphAntarctica*> grGrids; ///< The internally stored lat/lon grids.

    Int_t fGridPoints; ///< Number of points on the lat/lon grid.
    Bool_t needRemakeGrid; // if changed grid settings
    void updateHist();
    void init(RampdemReader::dataSet dataSet, Int_t coarseness); // in case I ever want another constructor

    void setPadMargins(); // prettification

  Bool_t fUseToolTip;
  TGToolTip* fToolTip;
  TString fToolTipUnits;
  void setToolTipUnits();

  Bool_t fBases;
  std::vector<TGraphAntarctica*> grBases;
  void updateBases();

  Bool_t fDrawnSelf; // Set by Draw(), help the updateHist() function to do sensible things
  void updateGPadPrims(std::vector<TGraphAntarctica*>& grs, Bool_t drawThem, Option_t* opt);
  ClassDef(AntarcticaBackground, 0)

  std::vector<Int_t> fOldPalette; //! Don't persist
  Bool_t fOldGrayScale; //! Don't persist
  TExec* fPalSetter; //! Don't persist
  TExec* fPalUnsetter; //! Don't persist
  std::vector<TPad*> fPads; //! Don't persist
};



#endif
