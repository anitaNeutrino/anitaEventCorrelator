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


ClassImp(AntarcticaBackground)



// simple global variable to increment integer in name so we can have as manny of these background as we like.
static int numAntarcticaBackgrounds = 0;


AntarcticaBackground::AntarcticaBackground(RampdemReader::dataSet dataSet, Int_t coarseness)
  : TProfile2D() {
  init(dataSet, coarseness);
}


AntarcticaBackground::~AntarcticaBackground(){
  deleteGrid();
}


void AntarcticaBackground::init(RampdemReader::dataSet dataSet, Int_t coarseness){
  SetDirectory(0);
  fName = Form("%s%d", getDefaultName(), numAntarcticaBackgrounds);
  numAntarcticaBackgrounds++;
  fDataSet = dataSet;
  fCoarseness = coarseness;
  needRemakeHist = true;  // force updateHist() to read in data by setting impossible coarseness
  fDrawnSelf = false; // set by Draw(), needed by updateHist()
  updateHist();

  fGridPoints = 1000;
  fDeltaLon = 30; // degrees
  fDeltaLat = 5; // degrees
  fGrid = true;
  needRemakeGrid = true; // force updateGrid() to make grid TGraphs on first call
  updateGrid();

  fUseToolTip = true;
  fToolTip = NULL;


  fExec = new TExec("fAntarcticaBackgroundExec",Form("%s->setPalette()", fName.Data()));


  // at some point, supporting ROOT versions < 6 is gonna be impossible...
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  palettes[RampdemReader::rampdem] = kLightTerrain;
  palettes[RampdemReader::bed] = kLake; // maybe a bit intense...
  palettes[RampdemReader::icemask_grounded_and_shelves] = kLightTerrain;
  palettes[RampdemReader::surface] = kLightTerrain;
  palettes[RampdemReader::thickness] = kRedBlue;
#endif
}



void AntarcticaBackground::setPalette(){

  // at some point, supporting ROOT versions < 6 is gonna be impossible...
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)

  // here I define the color palettes tbat I want the different data sets to use when they're drawn on the background.
  std::map<RampdemReader::dataSet, EColorPalette>::iterator it = palettes.find(fDataSet);
  if(it != palettes.end()){
    gStyle->SetPalette(it->second);
    gStyle->SetNdivisions(254);
  }
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

    if(fDrawnSelf){

      // now set the viewing range to the same as before if we have a gPad instance
      fXaxis.SetRangeUser(lowX, highX);
      fYaxis.SetRangeUser(lowY, highY);

      // and update the z-axis title if needed
      prettifyPalette();
    }

    setToolTipUnits();
  }

  // prettification
  GetXaxis()->SetNdivisions(0, kFALSE);
  GetYaxis()->SetNdivisions(0, kFALSE);

  needRemakeGrid = false;


}






Int_t AntarcticaBackground::GetCoarseness(){
  return fCoarseness;
}



void AntarcticaBackground::ToolTip(Bool_t toolTip){
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

Bool_t AntarcticaBackground::GetToolTip(){
  return fUseToolTip;
}


void AntarcticaBackground::SetCoarseness(Int_t coarseness){

  // sanity check
  if(coarseness < 1){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", coarsenesss must be >= 1. Setting coareness = "
	      << defaultCoarseness << std::endl;
    coarseness = defaultCoarseness;
  }

  needRemakeHist = fCoarseness == coarseness ? false : true;
  fCoarseness = coarseness;
  updateHist();
}


RampdemReader::dataSet AntarcticaBackground::GetDataSet(){
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

      // gr->SetDrawOption("*");
      gr->SetMarkerStyle(8);
      gr->SetMarkerColor(kMagenta);

      gr->SetDrawOption("p");
      // if(fBases){
      // 	gr->Draw("psame");
      // }

      grBases.push_back(gr);
    }
  }

  updateGPadPrims(grBases, fBases, "p");
}





void AntarcticaBackground::updateGrid(){

  if(needRemakeGrid){

    const int minLat = -90 + fDeltaLat;
    const int maxLat = -50;

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








void AntarcticaBackground::updateGPadPrims(std::vector<TGraphAntarctica*>& grs, Bool_t drawThem, Option_t* opt){

  if(gPad){
    TList* prims = gPad->GetListOfPrimitives();

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

      // now add the TGraphs with the proper options
      for(UInt_t grInd=0; grInd < grs.size(); grInd++){
	TGraphAntarctica* gr = grs.at(grInd);
	prims->AddLast(gr, opt);
      }

      // and re-add the things we removed
      for(UInt_t i=0; i < tempObjs.size(); i++){
	prims->AddLast(tempObjs.at(i), tempOpts.at(i));
      }
    }
    else{
      // Remove the TGraphsAntarcticas from the list of pad primitives
      for(UInt_t grInd=0; grInd < grs.size(); grInd++){
	TGraph* gr = grs.at(grInd);
	prims->RecursiveRemove(gr);
      }
    }
    gPad->Modified();
    gPad->Update();
  }
}





void AntarcticaBackground::deleteGrid(){

  // Grid(false); // avoid segfault?
  while(grGrids.size() > 0){

    TGraph* gr = grGrids.back();
    if(gPad){
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
 * @param opt is the draw option (default is colz)
 */
void AntarcticaBackground::Draw(Option_t* opt){

  TProfile2D::Draw(opt);
  fExec->Draw();
  TProfile2D::Draw(Form("%s same", opt));
  setPadMargins();
  prettifyPalette();
  fXaxis.SetAxisColor(kWhite);
  fYaxis.SetAxisColor(kWhite);

  gPad->Update();


  fDrawnSelf = true;

  SetBit(kMustCleanup);
  SetBit(kCanDelete); // This means the TPad that we've drawn on owns this, and should delete it when the TPad is destroyed.
  updateGrid();
  ToolTip(fToolTip);

}


/**
 * Helper function which prettifies the pad
 */
void AntarcticaBackground::setPadMargins(){
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.02);
  gPad->SetLeftMargin(0.02);
  gPad->SetRightMargin(0.1);
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
void AntarcticaBackground::prettifyPalette(){
  gPad->Modified();
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*) GetListOfFunctions()->FindObject("palette");
  if(palette){
    palette->SetX1NDC(0.91);
    palette->SetX2NDC(0.95);
    palette->SetY1NDC(0.55);
    palette->SetY2NDC(0.95);
    // palette->SetX1NDC(0.03);
    // palette->SetX2NDC(0.06);
    // palette->SetY1NDC(0.03);
    // palette->SetY2NDC(0.16);
    // palette->SetTitleSize(0.001);
    // palette->SetTitleOffset(0.1);

    TAxis* zAxis = GetZaxis();
    zAxis->SetTitle(RampdemReader::dataSetToAxisTitle(fDataSet));
    zAxis->SetTitleSize(0.001);
    zAxis->SetLabelSize(0.001);
    // std::cout << zAxis->GetTitleOffset() << std::endl;
    zAxis->SetTitleOffset(25);
    gPad->Modified();
    gPad->Update();
  }
}




// Interactive functions...

void AntarcticaBackground::Rampdem(bool useRampdem){
  SetDataSet(RampdemReader::rampdem);
}

Bool_t AntarcticaBackground::GetRampdem(){
  return fDataSet == RampdemReader::rampdem;
}

void AntarcticaBackground::Bed(bool useBed){
  SetDataSet(RampdemReader::bed);
};

Bool_t AntarcticaBackground::GetBed(){
  return fDataSet == RampdemReader::bed;
}

void AntarcticaBackground::Icemask(bool useIcemask){
  SetDataSet(RampdemReader::icemask_grounded_and_shelves);
}

Bool_t AntarcticaBackground::GetIcemask(){
  return fDataSet == RampdemReader::icemask_grounded_and_shelves;
}

void AntarcticaBackground::Surface(bool useSurface){
  SetDataSet(RampdemReader::surface);
}

Bool_t AntarcticaBackground::GetSurface(){
  return fDataSet == RampdemReader::surface;
}

void AntarcticaBackground::Thickness(bool useThickness){
  SetDataSet(RampdemReader::thickness);
}

Bool_t AntarcticaBackground::GetThickness(){
  return fDataSet == RampdemReader::thickness;
}

void AntarcticaBackground::Grid(Bool_t grid){
  fGrid = grid;

  updateGrid();
}

Bool_t AntarcticaBackground::GetGrid(){
  return fGrid;
}


void AntarcticaBackground::ShowBases(Bool_t showBases){
  fBases = showBases;
  updateBases();
}

Bool_t AntarcticaBackground::GetShowBases(){
  return fBases;
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







void AntarcticaBackground::updateToolTip(Int_t event, Int_t x, Int_t y, const char* extraInfo){

  if(fUseToolTip){
    Double_t easting = gPad->AbsPixeltoX(x);
    Double_t northing = gPad->AbsPixeltoY(y);
    Double_t val = GetBinContent(FindBin(easting, northing));
    // gPad->AbsPixeltoX(y);
    Double_t lon, lat;
    RampdemReader::EastingNorthingToLonLat(easting, northing, lon, lat, fDataSet);

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
