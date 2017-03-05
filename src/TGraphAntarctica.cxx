#include "TGraphAntarctica.h"
#include "TROOT.h"
#include "TVirtualPad.h"

ClassImp(TGraphAntarctica)


// Don't want this to be part of the class when we start writing things to files.
static AntarcticaBackground* fBackground = NULL;


void TGraphAntarctica::SetPoint(Int_t i, Double_t lon, Double_t lat){
  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  TGraph::SetPoint(i, easting, northing);
}








void TGraphAntarctica::Draw(Option_t* option){

  TString opt = option;
  opt.ToLower();

  // I'm going to go for defaulting to pretty-ish points?
  opt = opt.Length() == 0 ? "p" : opt;
  this->SetMarkerStyle(8);


  Bool_t drawAntarctica = false;

  if(!gPad){
    // No current canvas, so create one
    gROOT->MakeDefCanvas();
  }
  TList* prims = gPad->GetListOfPrimitives();

  // If "same" is specified... do not draw an AntarcticaBackground.
  if(opt.Contains("same")){
    opt.ReplaceAll("same", "");
    drawAntarctica = false;
  }

  // Else if "a" is specified but not in "same", then do draw AntarcticaBackground.
  else if(opt.Contains("a")){
    opt.ReplaceAll("a", "");
    drawAntarctica = true;
  }

  else{
    std::cout << "prims " << prims->GetEntries() << std::endl;
    if(prims->GetEntries() <= 1){
      drawAntarctica = true;
    }
  }

  if(drawAntarctica){
    fBackground = new AntarcticaBackground();
    fBackground->Draw();
  }
  else{
    for(int i=0; i < prims->GetEntries(); i++){
      TString primName = prims->At(i)->GetName();
      if(primName.Contains(AntarcticaBackground::getDefaultName())){
	fBackground = (AntarcticaBackground*) prims->At(i);
	break;
      }
    }
  }

  TGraph::Draw(opt);

}



void TGraphAntarctica::init(){
  doneConversion = false;
  convertArrays();
}









void TGraphAntarctica::convertArrays(){
  if(!doneConversion){
    for(Int_t i=0; i < GetN(); i++){
      RampdemReader::LonLatToEastingNorthing(fX[i], fY[i], fX[i], fY[i]);
    }
  }

  doneConversion = true;
}




/**
 * Interactive magic.
 *
 * @param event is the user interaction
 * @param x is the x-coordinate of the pixel under the mouse
 * @param y is the y-coordinate of the pixel under the mouse
 */
void TGraphAntarctica::ExecuteEvent(Int_t event, Int_t x, Int_t y){

  if(fBackground){
    fBackground->updateToolTip(event, x, y, "I'm a point!");
  }

  TGraph::ExecuteEvent(event, x, y);
}
