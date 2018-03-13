#include "TH2DAntarctica.h"
#include "RampdemReader.h"
#include "AntarcticaBackground.h"
#include "TGraphAntarctica.h"
#include "TRandom3.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TList.h"

ClassImp(TProfile2DAntarctica)
ClassImp(TH2DAntarctica)

const double palX1 = AntarcticaBackgroundDefaults::zAxisRightMargin;
const double palX2 = AntarcticaBackgroundDefaults::zAxisRightMargin + AntarcticaBackgroundDefaults::zAxisWidth;
const double palY1 = 1 - AntarcticaBackgroundDefaults::zAxisTopBottomMargin - AntarcticaBackgroundDefaults::zAxisHeight;
const double palY2 = 1 - AntarcticaBackgroundDefaults::zAxisTopBottomMargin;

TProfile2DAntarctica::TProfile2DAntarctica(Int_t nx, Int_t ny)
    : TProfile2D(), fAntarcticaBackground(NULL){

  AntarcticaBackground* b = getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));
  accept_stereographic = false;
}


TProfile2DAntarctica::TProfile2DAntarctica(const char* name, const char* title, Int_t nx, Int_t ny)
    : TProfile2D(), fAntarcticaBackground(NULL){

  SetNameTitle(name, title);
  AntarcticaBackground* b = getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));
  accept_stereographic = false;
}


TProfile2DAntarctica::TProfile2DAntarctica(const char* name, const char* title, const std::vector<Double_t>& x, const std::vector<Double_t>& y)
  : TProfile2D(), fAntarcticaBackground(NULL){

  SetNameTitle(name, title);
  if(x.size() > 1 && y.size() > 1){
    SetBins(x.size()-1, &x[0], y.size()-1, &y[0]);
  }
  else{
    AntarcticaBackground* b = getBackground();

    Int_t nx = b->GetNbinsX();
    Int_t ny = b->GetNbinsY();
    SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
	    ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));

  }
  accept_stereographic = false;
}


void TProfile2DAntarctica::Draw(Option_t* opt){

  TString opt2(opt);
  if(!opt2.Contains("same")){
    AntarcticaBackground* b = getBackground();
    b->Draw();
    b->SetBit(kCanDelete, false);
  }
  else{
    TList* l = gPad->GetListOfPrimitives();
    for(int i=0; i < l->GetEntries(); i++){
      TString objName = l->At(i)->GetName();
      if(objName.Contains("fAntarctica") && !objName.Contains("PalSetter")){
	if(fAntarcticaBackground && fAntarcticaBackground != (AntarcticaBackground*) l->At(i)){
	  delete fAntarcticaBackground;
	}
	fAntarcticaBackground = (AntarcticaBackground*) l->At(i);
	// break;
      }
    }
  }
  TString sameOpt = opt2.Contains("same") ? opt2 : TString::Format("%s same", opt);
  TProfile2D::Draw(sameOpt);
  ResetColorAxis();
  RescaleBackground();
  gPad->Modified();
}


Int_t TProfile2DAntarctica::Fill(Double_t lon, Double_t lat, Double_t val){

  Double_t easting, northing;
  if (accept_stereographic)
  {
    easting = lon;
    northing = lat;
  }
  else
  {
    RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  }
  return TProfile2D::Fill(easting, northing, val);
}


/** 
 * This is just to try out drawing one 2D histogram on top of another background
 * 
 * @param nTimes 
 */
void TProfile2DAntarctica::FillRandomly(Int_t nTimes){

  TRandom3 randy;

  for(int i=0; i < nTimes; i++){

    double lon = randy.Uniform(0, 360);
    double lat = randy.Uniform(-90, -60);

    Fill(lon, lat);
  }
}



/**
 * Helper function which prettifies the z-axis
 */
void TProfile2DAntarctica::ResetColorAxis(){
  gPad->Modified();
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*) GetListOfFunctions()->FindObject("palette");
  if(palette){
    palette->SetX1NDC(palX1);
    palette->SetX2NDC(palX2);
    palette->SetY1NDC(palY1);
    palette->SetY2NDC(palY2);
    
    TAxis* zAxis = GetZaxis();
    // zAxis->SetTitle();
    zAxis->SetTitleSize(AntarcticaBackgroundDefaults::zAxisTextSize);
    zAxis->SetLabelSize(AntarcticaBackgroundDefaults::zAxisTextSize);
    // std::cout << zAxis->GetTitleOffset() << std::endl;
    // zAxis->SetTitleOffset(0.1);
    gPad->Modified();
    gPad->Update();
  }
}


void TProfile2DAntarctica::ExecuteEvent(int event, int px, int py){
  if(event==kButton1Double){
    this->UnZoom();
  }
  AntarcticaBackground* b = getBackground();
  if(b){
    b->ExecuteEvent(event, px, py);
  }
  TProfile2D::ExecuteEvent(event, px, py);
}


TGraphAntarctica* TProfile2DAntarctica::findLocalMaxima() const {

  TGraphAntarctica* grLocalMaxima = new TGraphAntarctica();

  for(int by=2; by <= GetNbinsY()-1; by++){
    for(int bx=2; bx <= GetNbinsX()-1; bx++){

      // check the 8 points around this point
      double cVal = GetBinContent(bx, by);

      // Is this a local maximum?
      if(cVal > GetBinContent(bx-1, by-1) &&
         cVal > GetBinContent(bx-1, by  ) &&
         cVal > GetBinContent(bx-1, by+1) &&
         cVal > GetBinContent(bx  , by-1) &&
         cVal > GetBinContent(bx  , by+1) &&
         cVal > GetBinContent(bx+1, by-1) &&
         cVal > GetBinContent(bx+1, by  ) &&
         cVal > GetBinContent(bx+1, by+1)){

        Double_t easting  = fXaxis.GetBinCenter(bx);
        Double_t northing = fYaxis.GetBinCenter(by);
        grLocalMaxima->TGraph::SetPoint(grLocalMaxima->GetN(), easting, northing);
      }
    }
  }
  return grLocalMaxima;
}

TAxis* TProfile2DAntarctica::GetXaxis(){
  AntarcticaBackground* b = getBackground();
  if(b){
    return b->GetXaxis();
  }
  else{
    return &fXaxis;
  }
}

TAxis* TProfile2DAntarctica::GetYaxis(){
  AntarcticaBackground* b = getBackground();
  if(b){
    return b->GetYaxis();
  }
  else{
    return &fYaxis;
  }
}





















TH2DAntarctica::TH2DAntarctica(Int_t nx, Int_t ny)
    : TH2D(), fAntarcticaBackground(NULL){

  AntarcticaBackground* b = getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));

  accept_stereographic = false;
}


TH2DAntarctica::TH2DAntarctica(const char* name, const char* title, Int_t nx, Int_t ny)
    : TH2D(), fAntarcticaBackground(NULL){

  SetNameTitle(name, title);
  AntarcticaBackground* b = getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));
  accept_stereographic = false; 
}



TH2DAntarctica::TH2DAntarctica(const char* name, const char* title, const std::vector<Double_t>& x, const std::vector<Double_t>& y)
    : TH2D(), fAntarcticaBackground(NULL){

  SetNameTitle(name, title);
  if(x.size() > 1 && y.size() > 1){
    SetBins(x.size()-1, &x[0], y.size()-1, &y[0]);
  }
  else{
    AntarcticaBackground* b = getBackground();

    Int_t nx = b->GetNbinsX();
    Int_t ny = b->GetNbinsY();
    SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
	    ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));

  }
  accept_stereographic = false;
}




void TH2DAntarctica::Draw(Option_t* opt){

  TString opt2(opt);
  if(!opt2.Contains("same")){
    AntarcticaBackground* b = getBackground();
    b->Draw();
    b->SetBit(kCanDelete, false);
  }
  else{
    TList* l = gPad->GetListOfPrimitives();
    for(int i=0; i < l->GetEntries(); i++){
      TString objName = l->At(i)->GetName();
      if(objName.Contains("fAntarctica") && !objName.Contains("PalSetter")){
	if(fAntarcticaBackground && fAntarcticaBackground != (AntarcticaBackground*) l->At(i)){
	  delete fAntarcticaBackground;
	}
	fAntarcticaBackground = (AntarcticaBackground*) l->At(i);
	// break;
      }
    }
  }
  TString sameOpt = opt2.Contains("same") ? opt2 : TString::Format("%s same", opt);
  TH2D::Draw(sameOpt);
  ResetColorAxis();
  RescaleBackground();
  gPad->Modified();
}


Int_t TH2DAntarctica::Fill(Double_t lon, Double_t lat, Double_t val){

  Double_t easting, northing;

  if (accept_stereographic) 
  {
    easting = lon; 
    northing = lat; 
  }
  else
  {
    RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  }
  return TH2D::Fill(easting, northing, val);  
}


/** 
 * This is just to try out drawing one 2D histogram on top of another background
 * 
 * @param nTimes 
 */
void TH2DAntarctica::FillRandomly(Int_t nTimes){

  TRandom3 randy;

  for(int i=0; i < nTimes; i++){

    double lon = randy.Uniform(0, 360);
    double lat = randy.Uniform(-90, -60);

    Fill(lon, lat);
  }
}


/**
 * Helper function which prettifies the z-axis
 */
void TH2DAntarctica::ResetColorAxis(){
  if(gPad){
    gPad->Modified();
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis*) GetListOfFunctions()->FindObject("palette");
    if(palette){
      palette->SetX1NDC(palX1);
      palette->SetX2NDC(palX2);
      palette->SetY1NDC(palY1);
      palette->SetY2NDC(palY2);
      TAxis* zAxis = GetZaxis();
      // zAxis->SetTitle();
      zAxis->SetTitleSize(AntarcticaBackgroundDefaults::zAxisTextSize);
      zAxis->SetLabelSize(AntarcticaBackgroundDefaults::zAxisTextSize);
      // std::cout << zAxis->GetTitleOffset() << std::endl;
      // zAxis->SetTitleOffset(0.1);
      gPad->Modified();
      gPad->Update();
    }
  }
}


void TH2DAntarctica::ExecuteEvent(int event, int px, int py){
  if(event==kButton1Double){
    this->UnZoom();
  }
  AntarcticaBackground* b = getBackground();
  if(b){
    b->ExecuteEvent(event, px, py);
  }
  TH2D::ExecuteEvent(event, px, py);
}


TGraphAntarctica* TH2DAntarctica::findLocalMaxima() const {

  TGraphAntarctica* grLocalMaxima = new TGraphAntarctica();

  for(int by=2; by <= GetNbinsY()-1; by++){
    for(int bx=2; bx <= GetNbinsX()-1; bx++){

      // check the 8 points around this point
      double cVal = GetBinContent(bx, by);

      // Is this a local maximum?
      if(cVal > GetBinContent(bx-1, by-1) &&
         cVal > GetBinContent(bx-1, by  ) &&
         cVal > GetBinContent(bx-1, by+1) &&
         cVal > GetBinContent(bx  , by-1) &&
         cVal > GetBinContent(bx  , by+1) &&
         cVal > GetBinContent(bx+1, by-1) &&
         cVal > GetBinContent(bx+1, by  ) &&
         cVal > GetBinContent(bx+1, by+1)){

        Double_t easting  = fXaxis.GetBinCenter(bx);
        Double_t northing = fYaxis.GetBinCenter(by);
        grLocalMaxima->TGraph::SetPoint(grLocalMaxima->GetN(), easting, northing);
      }
    }
  }
  return grLocalMaxima;
}


TAxis* TH2DAntarctica::GetXaxis(){
  AntarcticaBackground* b = getBackground();
  if(b){
    return b->GetXaxis();
  }
  else{
    return &fXaxis;
  }
}

TAxis* TH2DAntarctica::GetYaxis(){
  AntarcticaBackground* b = getBackground();
  if(b){
    return b->GetYaxis();
  }
  else{
    return &fYaxis;
  }
}
