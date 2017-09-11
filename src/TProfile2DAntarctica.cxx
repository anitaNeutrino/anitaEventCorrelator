#include "TProfile2DAntarctica.h"
#include "RampdemReader.h"
#include "AntarcticaBackground.h"
#include "TRandom3.h"
ClassImp(TProfile2DAntarctica)


TProfile2DAntarctica::TProfile2DAntarctica(Int_t nx, Int_t ny)
: TProfile2D(), fAntarcticaBackground(NULL) {

  AntarcticaBackground* b =getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(nx),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(ny));
  
}


TProfile2DAntarctica::~TProfile2DAntarctica(){
  if(fAntarcticaBackground){
    delete fAntarcticaBackground;
    fAntarcticaBackground = NULL;
  }
}


AntarcticaBackground* TProfile2DAntarctica::getBackground(){
  if(!fAntarcticaBackground){
    fAntarcticaBackground = new AntarcticaBackground();
  }
  return fAntarcticaBackground;
}


void TProfile2DAntarctica::Draw(Option_t* opt){

  AntarcticaBackground* b = getBackground();
  b->SetCoarseness(25);
  b->SetNdivisions(10);
  b->Draw("colz");
  // should have constructed a new TPad, so can omit same...
  
  TString sameOpt = TString::Format("%s same", opt);
  TProfile2D::Draw(sameOpt); 
}



Int_t TProfile2DAntarctica::Fill(Double_t lon, Double_t lat, Double_t weight){

  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  return TProfile2D::Fill(easting, northing, weight);  
}


/** 
 * This is just to try out drawing one 2D histogram on top of another background
 * 
 * @param nTimes 
 */
void TProfile2DAntarctica::FillRandom(Int_t nTimes){

  TRandom3 randy;

  for(int i=0; i < nTimes; i++){

    double lon = randy.Uniform(0, 360);
    double lat = randy.Uniform(-90, -60);

    Fill(lon, lat);
  }
  
  
  
}
