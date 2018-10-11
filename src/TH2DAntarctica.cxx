#include "TH2DAntarctica.h"
#include "RampdemReader.h"
#include "AntarcticaBackground.h"
#include "TGraphAntarctica.h"
#include "TRandom3.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TList.h"
#include "AntarcticaGeometry.h"

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

void getAngularResolution(Double_t snr, Double_t& sigmaTheta, Double_t& sigmaPhi)
{
  //this is just here for the error contour filling
  sigmaPhi = (AnitaVersion::get() == 3) ? exp(-2.50414e-01 * snr + 3.02406e-01) + 2.43376e-01 : 5.09057/(pow(snr, 8.01369e-01) + 1.);
  sigmaTheta = (AnitaVersion::get() == 3) ? exp(-3.83773e-01 * snr + -3.00964e-01) + 1.64537e-01 : 1.34307/(pow(snr, 7.09382e-01) + 1.);
}

Int_t TH2DAntarctica::FillWithErrorContours(Double_t lon, Double_t lat, Double_t phi, Double_t theta,Double_t snr, Double_t ll_thresh, UsefulAdu5Pat upat, Double_t dist_thresh){
  Double_t sigmaPhi, sigmaTheta;
  getAngularResolution(snr, sigmaTheta, sigmaPhi);
  return FillWithErrorContours(lon,lat,phi,theta,sigmaPhi,sigmaTheta,ll_thresh,upat,dist_thresh);
}


Int_t TH2DAntarctica::FillWithErrorContours(Double_t lon, Double_t lat, Double_t phi, Double_t theta,Double_t sigmaPhi, Double_t sigmaTheta, Double_t ll_thresh, UsefulAdu5Pat upat, Double_t dist_thresh){
  //maxval is the value of the reconstructed point on the continent, other points are maxval - sqrt(ll)
  Int_t maxval = ceil(sqrt(ll_thresh)) + 1;

  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);

  Int_t binX = TH2D::GetXaxis()->FindBin(easting);
  Int_t binY = TH2D::GetYaxis()->FindBin(northing);

  bool checking = true;
  Int_t n_added = 1;

  //want to keep track of which ways it's still worth checking
  int udlr[4]; 
  bool updone = 0;
  bool downdone = 0;
  bool leftdone = 0;
  bool rightdone = 0;
  Int_t nrow;
  Int_t ncol;
  int level = 1;

  while(checking)
  {
    int miss_count = 0;
    int curr = 0;
    udlr[0] = 0;
    udlr[1] = 0;
    udlr[2] = 0;
    udlr[3] = 0;
    for(int i = 0; i <  8*level; i++)
    {
      Double_t thetaSource, phiSource;
      Double_t sourceEasting, sourceNorthing;
      Double_t sourceLon, sourceLat;
      nrow = 3 + 2*(level-1);
      ncol = nrow - 2;

      Int_t newBinX, newBinY;

      //search outwards from a point in squares
      if(i < nrow){
        if(updone){
          miss_count++;
          continue;
        }
        curr = 0;
        newBinX = binX-level+i;
        newBinY = binY+level;
      }
      else if(i < 2 * nrow){
        if(downdone){
          miss_count++;
          continue;
        }
        curr = 1;
        newBinX = binX-level+i-nrow;
        newBinY = binY-level;
      }
      else{
        Int_t alternator = i%2 ? 1 : -1;
        if(alternator == 1 && rightdone){
          miss_count++;
          continue;
        }
        if(alternator == -1 && leftdone){
          miss_count++;
          continue;
        }
        curr = alternator == 1 ? 3 : 2;
        Int_t ypos = level - 1 - (i - 2*nrow)/2;
        newBinX = binX+(alternator*level);
        newBinY = binY+ypos;
      }
      sourceEasting = TH2D::GetXaxis()->GetBinCenter(newBinX);
      sourceNorthing = TH2D::GetYaxis()->GetBinCenter(newBinY);
      RampdemReader::EastingNorthingToLonLat(sourceEasting, sourceNorthing, sourceLon, sourceLat);
      Double_t sourceAlt = RampdemReader::SurfaceAboveGeoid(sourceLon, sourceLat);

      if(dist_thresh > 0)
      {
        Double_t surfaceDist = 1e-3*AntarcticCoord::surfaceDistance(sourceLat, lat, sourceLon, lon);
        if(surfaceDist < dist_thresh)
        {
          Fill(sourceLon, sourceLat, maxval-1);
          n_added++;
          continue;
        }
      }
      
      upat.getThetaAndPhiWave2(sourceLon, sourceLat, sourceAlt, thetaSource, phiSource);
      
      thetaSource *= TMath::RadToDeg();
      phiSource *= TMath::RadToDeg();
      
      Double_t dTheta = (thetaSource - theta)/sigmaTheta;
      Double_t deltaPhi = phiSource - phi;
      if(deltaPhi > 180) {
        while(deltaPhi >= 180)
          deltaPhi -= 360;
        while(deltaPhi < -180)
          deltaPhi += 360;
      }
      Double_t dPhi = deltaPhi/sigmaPhi;

      Double_t ll = dTheta * dTheta + dPhi * dPhi;
      if(ll > ll_thresh)
      {
        miss_count++;
        udlr[curr]++;
        continue;
      }
      Fill(sourceLon, sourceLat, maxval - ceil(sqrt(ll)));
      n_added++;
    }
    if(miss_count == 8*level) checking = false;
    if(udlr[0] == nrow) updone = true;
    if(udlr[1] == nrow) downdone = true;
    if(udlr[2] == ncol) leftdone = true;
    if(udlr[3] == ncol) rightdone = true;
    level++;
  }
  Fill(lon, lat, maxval);
  return n_added;
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
