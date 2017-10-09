#include "AntarcticaBackground.h"
#include "TGraphAntarctica.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "BaseList.h"
#include "TStyle.h"
#include "TExec.h"
#include "TGToolTip.h"
#include "TColor.h"


ClassImp(AntarcticaBackground)



// simple global variable to increment integer in name so we can have as manny of these background as we like.
static int numAntarcticaBackgrounds = 0;


AntarcticaBackground::AntarcticaBackground(RampdemReader::dataSet dataSet, Int_t coarseness)
    : TProfile2D(), hDummy() {
  init(dataSet, coarseness);
}


AntarcticaBackground::~AntarcticaBackground(){
  if(gPad){
    SetToolTip(false);
    gPad->Modified();
    gPad->Update();
  }
  deleteGrid();
}


void AntarcticaBackground::init(RampdemReader::dataSet dataSet, Int_t coarseness){
  SetDirectory(0);
  fName = Form("%s%d", getDefaultName(), numAntarcticaBackgrounds);
  numAntarcticaBackgrounds++;
  fMinVal = DBL_MAX;
  fMaxVal = -DBL_MAX;
  fDataSet = dataSet;
  fCoarseness = coarseness;
  needRemakeHist = true;  // force updateHist() to read in data by setting impossible coarseness
  fDrawnSelf = false; // set by Draw(), needed by updateHist()
  updateHist();

  fGrayScale = false;  

  fGridPoints = 1000;
  fDeltaLon = 30; // degrees
  fDeltaLat = 5; // degrees
  fGrid = true;
  needRemakeGrid = true; // force updateGrid() to make grid TGraphs on first call
  updateGrid();

  fUseToolTip = true;
  fToolTip = NULL;

  TString palSetterName = TString::Format("%sPalSetter", fName.Data());
  fPalSetter = new TExec(palSetterName,Form("%s->setPalette()", fName.Data()));
  TString palUnsetterName = TString::Format("%sPalSetter", fName.Data());  
  fPalUnsetter = new TExec(palUnsetterName,Form("%s->unsetPalette()", fName.Data()));

  // Dummy histogram whose only purpose is to have a color axis, which we will steal
  // and stay well away from out Antarctica histograms, hence the crazy axis limits
  hDummy.SetBins(2, -9e30, -8e30, 1, -9e30, -8e30);
  fShowColorAxis = true;

  // at some point, supporting ROOT versions < 6 is gonna be impossible...
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  palettes[RampdemReader::rampdem] = kLightTerrain;
  palettes[RampdemReader::bed] = kLake; // maybe a bit intense...
  palettes[RampdemReader::icemask_grounded_and_shelves] = kLightTerrain;
  palettes[RampdemReader::surface] = kLightTerrain;
  palettes[RampdemReader::thickness] = kRedBlue;
  opacity = 1; 
#endif

  fBases = false;
  fDrawBasesOnTop = true;
}



void AntarcticaBackground::setPalette(){

  // at some point, supporting ROOT versions < 6 is gonna be impossible...
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)

  // here I define the color palettes tbat I want the different data sets to use when they're drawn on the background.
  std::map<RampdemReader::dataSet, EColorPalette>::iterator it = palettes.find(fDataSet);

  if(it != palettes.end()){

    // save the old palette
    fOldPalette.resize(gStyle->GetNumberOfColors());
    for(UInt_t i=0; i < fOldPalette.size(); i++){
      fOldPalette[i] = gStyle->GetColorPalette(i);
    }
    fOldGrayScale = TColor::IsGrayscale();
    
    TColor::SetGrayscale(fGrayScale);
    gStyle->SetPalette(it->second,0,opacity);
  }
#else
  std::cerr << __PRETTY_FUNCTION__ << " requires ROOT version at least 6" << std::endl;
#endif

}


void AntarcticaBackground::unsetPalette(){

  // at some point, supporting ROOT versions < 6 is gonna be impossible...
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  TColor::SetGrayscale(fOldGrayScale);
  if(fOldPalette.size() > 0){
    gStyle->SetPalette(fOldPalette.size(), &fOldPalette[0]);
  }
  // std::cerr << "here" << std::endl;
    
#else
  std::cerr << __PRETTY_FUNCTION__ << " requires ROOT version at least 6" << std::endl;
#endif
}




/**
 * Update the map of Antarctica as different options are selected
 */
void AntarcticaBackground::updateHist(){


  if(needRemakeHist){

    // Here I save the current drawn axis range so I can set it to be the same after updating contents
    Int_t firstX = fXaxis.GetFirst();
    Int_t lastX = fXaxis.GetLast();
    Int_t firstY = fYaxis.GetFirst();
    Int_t lastY = fYaxis.GetLast();

    Double_t lowX = fXaxis.GetBinLowEdge(firstX);
    Double_t highX = fXaxis.GetBinUpEdge(lastX);
    Double_t lowY = fXaxis.GetBinLowEdge(firstY);
    Double_t highY = fXaxis.GetBinUpEdge(lastY);

    // std::cout << firstX << "\t" << lastX << "\t" << lowX << "\t" << highX << std::endl;

    // now I get the new histogram binning
    Int_t nx, ny;
    RampdemReader::getNumXY(nx, ny, fDataSet);
    Double_t xMin, xMax, yMin, yMax;
    RampdemReader::getMapCoordinates(xMin, yMin, xMax, yMax, fDataSet);

    // accounting for coarseness
    nx /= fCoarseness;
    ny /= fCoarseness;


    // Get rid of old the AntarcticaBackground content
    fBinEntries.Reset();
    for(int by=0; by <= GetNbinsY() + 1; by++){
      for(int bx=0; bx <= GetNbinsX() + 1; bx++){
        SetBinContent(bx, by, 0);
      }
    }
    // change the histogram dimensions
    SetBins(nx, xMin, xMax, ny, yMin, yMax);

    // insert new data
    RampdemReader::fillThisHist(this, fDataSet);

    fMinVal = DBL_MAX;
    fMaxVal = -DBL_MAX;
    for(int by=1; by <= GetNbinsY(); by++){
      for(int bx=1; bx <= GetNbinsX(); bx++){
        double val = GetBinContent(bx, by);
        if(val < fMinVal){
          fMinVal = val;
        }
        if(val > fMaxVal){
          fMaxVal = val;
        }
      }
    }
    SetMaximum(fMaxVal);
    SetMinimum(fMinVal);
    fScaleMax = 0;
    fScaleMin = 0;
    hDummy.SetBinContent(1, fMinVal);
    hDummy.SetBinContent(2, fMinVal);
    hDummy.SetMaximum(fMaxVal);
    hDummy.SetMinimum(fMinVal);

    if(fDrawnSelf){

      // now set the viewing range to the same as before if we have a gPad instance
      fXaxis.SetRangeUser(lowX, highX);
      fYaxis.SetRangeUser(lowY, highY);

      // and update the z-axis title if needed
      ResetColorAxis();
    }

    setToolTipUnits();
  }

  // prettification
  GetXaxis()->SetNdivisions(0, kFALSE);
  GetYaxis()->SetNdivisions(0, kFALSE);

  needRemakeGrid = false;
}






Int_t AntarcticaBackground::GetCoarseness() const{
  return fCoarseness;
}



void AntarcticaBackground::SetToolTip(Bool_t toolTip){
  fUseToolTip = toolTip;

  if(fToolTip && !fUseToolTip){
    delete fToolTip;
    fToolTip = NULL;
  }
  else{
    if(!fToolTip){
      fToolTip = new TGToolTip();
      fToolTip->SetBackgroundColor(kWhite);
    }
  }
}

Bool_t AntarcticaBackground::GetToolTip() const{
  return fUseToolTip;
}


void AntarcticaBackground::SetCoarseness(Int_t coarseness){

  // sanity check
  if(coarseness < 1){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", coarsenesss must be >= 1. Setting coareness = "
              << AntarcticaBackgroundDefaults::defaultCoarseness << std::endl;
    coarseness = AntarcticaBackgroundDefaults::defaultCoarseness;
  }

  needRemakeHist = fCoarseness == coarseness ? false : true;
  fCoarseness = coarseness;
  updateHist();
}


RampdemReader::dataSet AntarcticaBackground::GetDataSet() const{
  return fDataSet;
}



void AntarcticaBackground::SetDataSet(RampdemReader::dataSet dataSet){

  needRemakeHist = dataSet == fDataSet ? false : true;
  fDataSet = dataSet;
  updateHist();
}







void AntarcticaBackground::updateBases(){

  if(grBases.size()==0){
    for(UInt_t b=0; b < BaseList::getNumBases(); b++){
      TGraphAntarctica* gr = new TGraphAntarctica();
      BaseList::base theBase = BaseList::getBase(b);
      gr->SetPoint(gr->GetN(), theBase.longitude, theBase.latitude);
      gr->SetTitle(theBase.name);
      gr->SetName(theBase.name);

      // gr->SetDrawOption("*");
      gr->SetMarkerStyle(8);
      gr->SetMarkerColor(kMagenta);
      gr->SetEditable(false);

      gr->SetDrawOption("p");
      // if(fBases){
      // 	gr->Draw("psame");
      // }

      grBases.push_back(gr);
    }
  }

  updateGPadPrims(grBases, fBases, "p", fDrawBasesOnTop);
}


void AntarcticaBackground::updateGrid(){

  if(needRemakeGrid){

    const int minLat = -90 + fDeltaLat;
    const int maxLat = -60;

    deleteGrid();

    // make circles of constant latitude
    for(Int_t lat = minLat; lat<= maxLat; lat += fDeltaLat){
      TGraphAntarctica* gr = new TGraphAntarctica();
      gr->SetLineColor(kGray);
      const Double_t deltaLon = 360./fGridPoints;
      for(int i=0; i < fGridPoints; i++){
        Double_t theLat = lat;
        Double_t theLon = i*deltaLon;
        gr->SetPoint(gr->GetN(), theLon, theLat);
        // std::cout << gr << "\t" << gr->GetN() << "\t" << easting << "\t" << northing << std::endl;
      }
      gr->SetEditable(false);
      gr->SetTitle(Form("Grid: Lat %d", lat)); // descriptive title
      grGrids.push_back(gr);
    }

    // make lines of constant longitude
    for(Int_t lon = -180; lon < 180; lon+= fDeltaLon){
      TGraphAntarctica* gr = new TGraphAntarctica();
      gr->SetLineColor(kGray);
      const Double_t deltaLat = double(maxLat - -90)/fGridPoints;
      for(int i=0; i < fGridPoints; i++){
        Double_t theLat = -90 + deltaLat*i;
        Double_t theLon = lon;
        gr->SetPoint(gr->GetN(), theLon, theLat);
      }
      gr->SetEditable(false);
      gr->SetTitle(Form("Grid: Lon %d", lon)); // descriptive title
      grGrids.push_back(gr);
    }

    needRemakeGrid = false;
  }

  updateGPadPrims(grGrids, fGrid, "l");
}








void AntarcticaBackground::updateGPadPrims(std::vector<TGraphAntarctica*>& grs, Bool_t drawThem, Option_t* opt, bool drawGraphsOnTop){


  // for(UInt_t padInd = 0; padInd < fPads.size(); padInd++){
  TVirtualPad* fPad = gPad;
  if(fPad){
    TList* prims = fPad->GetListOfPrimitives();

    if(drawThem && fDrawnSelf){
      // Manually add the grid TGraphs into the list of pad primitives...
      // turns out that to do this with options is slightly more complicated than I thought.
      // You need to iterate over the links (that wrap the objects) because the links hold the options.
      TObjLink *thisLink = prims->FirstLink();

      Int_t numThis = 0;
      // while(thisLink->GetObject() != this){
      while(numThis < 2){
        thisLink = thisLink->Next();

        // Complication, now we're drawing two copies of the histogram with an exec to set the palette
        // I need to find the second instance on the canvas
        if(thisLink->GetObject()==this){
          numThis++;
        }
      }

      // so thisLink points to this.
      // using that info with add After.
      TObjLink* antarcticaStuff = thisLink->Next();

      // It turns out that for such a fundamental ROOT container, TList is pretty dumb.
      // I think I'm going to have to take all the objects out, save their associated options,
      // and redraw them. That's becasuse there's no way to add a link with options at an arbitrary
      // position in the list, only at the end of the list.

      // first copy the pointers...
      std::vector<TObject*> tempObjs;
      std::vector<TString> tempOpts;
      while(antarcticaStuff){
        tempObjs.push_back(antarcticaStuff->GetObject());
        tempOpts.push_back(antarcticaStuff->GetOption());
        antarcticaStuff = antarcticaStuff->Next();
      }

      // now take the objects out
      for(UInt_t i=0; i < tempObjs.size(); i++){
        prims->RecursiveRemove(tempObjs.at(i));
      }


      // now we need to redraw them in the proper order specified by onTop
      
      if(drawGraphsOnTop){
        // first, re-add the things we removed
        for(UInt_t i=0; i < tempObjs.size(); i++){
          prims->AddLast(tempObjs.at(i), tempOpts.at(i));
        }
        // then the TGraphs, with the proper options
        for(UInt_t grInd=0; grInd < grs.size(); grInd++){
          TGraphAntarctica* gr = grs.at(grInd);
          prims->AddLast(gr, opt);
        }
      }
      else{
        // otherwise add the TGraphs first, with the proper options
        for(UInt_t grInd=0; grInd < grs.size(); grInd++){
          TGraphAntarctica* gr = grs.at(grInd);
          prims->AddLast(gr, opt);
        }
        // and then re-add the things we removed
        for(UInt_t i=0; i < tempObjs.size(); i++){
          prims->AddLast(tempObjs.at(i), tempOpts.at(i));
        }
      }
    }
    else{
      // Remove the TGraphsAntarcticas from the list of pad primitives
      for(UInt_t grInd=0; grInd < grs.size(); grInd++){
        TGraph* gr = grs.at(grInd);
        prims->RecursiveRemove(gr);
      }
    }
    fPad->Modified();
    fPad->Update();
  }
  // }
}




void AntarcticaBackground::deleteGrid(){

  // Grid(false); // avoid segfault?
  while(grGrids.size() > 0){

    TGraph* gr = grGrids.back();
    // for(UInt_t padInd=0; padInd < fPads.size(); padInd++){
    TVirtualPad* fPad = gPad; //fPads[padInd];
    if(fPad){
      TList* prims = gPad->GetListOfPrimitives();
      // prims->Remove(gr);
      prims->RecursiveRemove(gr);
    }

    delete gr;
    grGrids.pop_back();
  }

}




void AntarcticaBackground::SetGridDivisions(Int_t deltaLon, Int_t deltaLat){

  needRemakeGrid = fDeltaLon == deltaLon && fDeltaLat == deltaLat ? false : true;
  fDeltaLon = deltaLon;
  fDeltaLat = deltaLat;
  updateGrid();
}



/**
 * Draw function for Antarctica Background which also prettifies the pad/canvas.
 *
 * @param opt is the draw option (default is col)
 */
void AntarcticaBackground::Draw(Option_t* opt){

  if(!gPad){
    gROOT->MakeDefCanvas();
  }

  // handle defaults
  TString opt2 = opt;
  if(opt2 == AntarcticaBackgroundDefaults::drawOpt){
    opt2 = "colz";
    fShowColorAxis = true;
  }
  if(opt2.Contains("z")){
    fShowColorAxis = true;
  }
  else{
    fShowColorAxis = false;
  }

  // under no circumstances allow this to draw it's own color axis...
  TString opt3 = opt2;
  opt3.ReplaceAll("z", "");

  // 1st instance, no color axis
  TProfile2D::Draw(opt3);

  fPalSetter->Draw(); // 1st exec, to set the background palette

  // Now use hDummy's palette instead, which won't be affected by rescaling
  if(fShowColorAxis){
    hDummy.Draw("colz same");
  }
  else{
    hDummy.Draw("col");
  }

  // 2nd instance, no color axis
  opt3 += " same";
  TProfile2D::Draw(opt3);

  // 2nd exec, to set the foreground palette again
  fPalUnsetter->Draw();


  // pad prettification
  setPadMargins();
  ResetColorAxis(true);
  fXaxis.SetAxisColor(kWhite);
  fYaxis.SetAxisColor(kWhite);
  fXaxis.SetTitleOffset(9999);
  fYaxis.SetTitleOffset(9999);  
  fXaxis.SetLabelOffset(9999);
  fYaxis.SetLabelOffset(9999);    
 

  // force redraw
  gPad->Update();
  fDrawnSelf = true;

  SetBit(kMustCleanup);
  SetBit(kCanDelete); // This means the TPad that we've drawn on owns this, and should delete it when the TPad is destroyed.
  // (These bits are UNSET in TH2DAntarctica, which can own it's own background)

  // set up things the background owns, if requested
  updateGrid();
  SetToolTip(fToolTip);

}


/**
 * Helper function which prettifies the pad
 */
void AntarcticaBackground::setPadMargins(){
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.02);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.02);
  // gPad->SetTopMargin(0.05);
  // gPad->SetBottomMargin(0.05);
  // gPad->SetLeftMargin(0.05);
  // gPad->SetRightMargin(0.1);
  // gPad->SetRightMargin(0.05);
  gPad->SetFrameLineColor(0);
  gPad->SetFrameLineWidth(0);
  gPad->SetFrameBorderSize(0);

}



/**
 * Helper function which prettifies the z-axis
 */
void AntarcticaBackground::ResetColorAxis(bool trigger_redraw){

  if(fShowColorAxis){
    if(trigger_redraw){
      gPad->Modified();
      gPad->Update();
    }

    // now, use the Dummy histogram's color axis
    TPaletteAxis *palette = (TPaletteAxis*) hDummy.GetListOfFunctions()->FindObject("palette");
    if(palette){
      palette->SetX1NDC(AntarcticaBackgroundDefaults::zAxisRightMargin);
      palette->SetX2NDC(AntarcticaBackgroundDefaults::zAxisRightMargin + AntarcticaBackgroundDefaults::zAxisWidth);
      palette->SetY1NDC(AntarcticaBackgroundDefaults::zAxisTopBottomMargin);
      palette->SetY2NDC(AntarcticaBackgroundDefaults::zAxisTopBottomMargin + AntarcticaBackgroundDefaults::zAxisHeight);

      TAxis* zAxis = GetZaxis();
      zAxis->SetTitle(RampdemReader::dataSetToAxisTitle(fDataSet));
      zAxis->SetTitleSize(AntarcticaBackgroundDefaults::zAxisTextSize);
      zAxis->SetLabelSize(AntarcticaBackgroundDefaults::zAxisTextSize);

      TAxis* zAxis2 = hDummy.GetZaxis();
      zAxis2->SetTitle(RampdemReader::dataSetToAxisTitle(fDataSet));
      zAxis2->SetTitleSize(AntarcticaBackgroundDefaults::zAxisTextSize);
      zAxis2->SetLabelSize(AntarcticaBackgroundDefaults::zAxisTextSize);

      // std::cout << zAxis->GetTitleOffset() << std::endl;
      // zAxis->SetTitleOffset(0.1);
      if(trigger_redraw){
        gPad->Modified();
        gPad->Update();
      }
    }
    else{
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find dummy color axis! " << std::endl;
    }
  }
}



// Interactive functions...

void AntarcticaBackground::SetRampdem(bool useRampdem){
  SetDataSet(RampdemReader::rampdem);
}

Bool_t AntarcticaBackground::GetRampdem() const {
  return fDataSet == RampdemReader::rampdem;
}

void AntarcticaBackground::SetBed(bool useBed){
  SetDataSet(RampdemReader::bed);
};

Bool_t AntarcticaBackground::GetBed() const {
  return fDataSet == RampdemReader::bed;
}

void AntarcticaBackground::SetIcemask(bool useIcemask){
  SetDataSet(RampdemReader::icemask_grounded_and_shelves);
}

Bool_t AntarcticaBackground::GetIcemask() const {
  return fDataSet == RampdemReader::icemask_grounded_and_shelves;
}

void AntarcticaBackground::SetSurface(bool useSurface){
  SetDataSet(RampdemReader::surface);
}

Bool_t AntarcticaBackground::GetSurface() const {
  return fDataSet == RampdemReader::surface;
}

void AntarcticaBackground::SetThickness(bool useThickness){
  SetDataSet(RampdemReader::thickness);
}

Bool_t AntarcticaBackground::GetThickness() const {
  return fDataSet == RampdemReader::thickness;
}

void AntarcticaBackground::SetGrid(Bool_t grid){
  fGrid = grid;

  updateGrid();
}

Bool_t AntarcticaBackground::GetGrid() const {
  return fGrid;
}


void AntarcticaBackground::SetShowBases(Bool_t showBases){
  fBases = showBases;
  updateBases();
}

Bool_t AntarcticaBackground::GetShowBases() const {
  return fBases;
}


void AntarcticaBackground::SetGrayScale(bool grayScale){
  fGrayScale = grayScale;
}

Bool_t AntarcticaBackground::GetGrayScale() const {
  return fGrayScale;
}


void AntarcticaBackground::SetShowColorAxis(bool showColorAxis){
  fShowColorAxis = showColorAxis;
  if(fShowColorAxis){
    hDummy.SetDrawOption("colz same");
    ResetColorAxis(true);
  }
  else{
    hDummy.SetDrawOption("col same");
  }
}



/**
 * Interactive magic.
 *
 * @param event is the user interaction
 * @param x is the x-coordinate of the pixel under the mouse
 * @param y is the y-coordinate of the pixel under the mouse
 */
void AntarcticaBackground::ExecuteEvent(Int_t event, Int_t x, Int_t y){

  updateToolTip(event, x, y);
  TProfile2D::ExecuteEvent(event, x, y);
}



/** 
 * This is the magic function that makes two 2D histograms plotted on top of one 
 * another have the same dynamic range, by rescaling the background.
 * This is called by TH2DAntarctica using G
 * 
 * @param newMin should be TH2DAntarctica::GetMinimum()
 * @param newMax should be TH2DAntarctica::GetMaximum()
 */
void AntarcticaBackground::scale(double newMin, double newMax){

  // This line prevents the empty bins in the TH2D filling the whole histogram
  if(newMin == 0){
    // epsilon
    newMin += 1e-15;
  }

  double currentMax = -DBL_MAX;
  double currentMin = DBL_MAX;

  if(fScaleMax == fScaleMin){
    // need to discover the current scale max/min
    // not sure we can get here now?

    for(int by=1; by <= GetNbinsY(); by++){
      for(int bx=1; bx <= GetNbinsX(); bx++){

        Int_t bin = GetBin(bx, by);
        double binEntries = GetBinEntries(bin);
        // cleverer way of checking this was a filled background bin
        if(binEntries > 0){
          double val = GetBinContent(bx, by);

          currentMin = val < currentMin ? val : currentMin;
          currentMax = val > currentMax ? val : currentMax;
        }
      }
    }
  }
  else{
    currentMax = fScaleMax;
    currentMin = fScaleMin;
  }

  double currentDelta = currentMax - currentMin;
  double newDelta = newMax - newMin;

  for(int by=1; by <= GetNbinsY(); by++){
    for(int bx=1; bx <= GetNbinsX(); bx++){

      // cleverer way of checking this was a filled background bin
      Int_t bin = GetBin(bx, by);
      double binEntries = GetBinEntries(bin);

      if(binEntries > 0){

        double val = GetBinContent(bx, by);
        double frac = (val - currentMin)/currentDelta;

        // if(val != 0 && currentMax - val < 100){
        //   std::cerr << val << "\t" << currentMin << "\t" << currentMax << "\t" << frac << std::endl;
        // }

        double newVal = newMin + frac*newDelta;
        SetBinContent(bx, by, newVal);
        SetBinEntries(bin, 1); // otherwise TProfile scales the value

      }
    }
  }
  // record new scale
  fScaleMax = newMax;
  fScaleMin = newMin;

  // set max/min
  SetMaximum(newMax);
  SetMinimum(newMin);

  // trigger re-paint of canvas
  gPad->Modified();
  gPad->Update();

}





void AntarcticaBackground::updateToolTip(Int_t event, Int_t x, Int_t y, const char* extraInfo){

  if(fUseToolTip){
    Double_t easting = gPad->AbsPixeltoX(x);
    Double_t northing = gPad->AbsPixeltoY(y);
    Double_t val = GetBinContent(FindBin(easting, northing));
    // gPad->AbsPixeltoX(y);
    Double_t lon, lat;
    RampdemReader::EastingNorthingToLonLat(easting, northing, lon, lat);

    TString theToolTipText = Form("Lon %4.2lf\nLat %4.2lf\n%4.2f%s", lon, lat, val, fToolTipUnits.Data());

    if(extraInfo!=NULL){
      theToolTipText += TString::Format("\n%s", extraInfo);
    }

    fToolTip->SetText(theToolTipText.Data());
    TCanvas* theCan = gPad->GetCanvas();
    Int_t topX = theCan->GetWindowTopX();
    Int_t topY = theCan->GetWindowTopY();

    const int xOffset = 10;
    const int yOffset = 10 + fToolTip->GetHeight()/2;


    fToolTip->Show(topX + x + xOffset, topY + y + yOffset);

  }
  TProfile2D::ExecuteEvent(event, x, y);
}




/**
 * Parses the z-axis title and extracts the units so they can be put in the tool tip.
 */
void AntarcticaBackground::setToolTipUnits(){
  TString tempAxisTitle(RampdemReader::dataSetToAxisTitle(fDataSet));
  TObjArray* tokens = tempAxisTitle.Tokenize("(");

  Bool_t gotUnits = false;
  int nTokens = tokens->GetEntries();
  // std::cout << tempAxisTitle << "\t" << nTokens << std::endl;
  if(nTokens > 1){
    TString afterOpenParen = ((TObjString*) tokens->At(1))->GetString();
    TObjArray* tokens2 = afterOpenParen.Tokenize(")");
    int nTokens2 = tokens2->GetEntries();
    // std::cout << afterOpenParen << "\t" << nTokens2 << std::endl;
    if(nTokens2 > 0){
      fToolTipUnits = " (" + TString(((TObjString*) tokens2->At(0))->GetString()) + ")";
      gotUnits = true;
    }
    // std::cout << "I set it to be " << fToolTipUnits.Data() << std::endl;
    delete tokens2;
  }
  delete tokens;

  if(!gotUnits){
    fToolTipUnits = "";
  }
}
