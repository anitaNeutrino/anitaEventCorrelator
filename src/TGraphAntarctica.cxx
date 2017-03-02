#include "TGraphAntarctica.h"
// #include "TROOT.h"
#include "TVirtualPad.h"



void TGraphAntarctica::Draw(Option_t* option){


  TString opt = option;

  // set the default drawing option to match TGraph
  opt = opt.Length() > 0 ? opt : "alp";
  opt.ToLower();

  // same should override a

  Bool_t drawAntarctica = false;

  if(opt.Contains("same")){
    opt.ReplaceAll("same", "");
  }
  else if(opt.Contains("a")){
    // then redraw background hist if needed
    opt.ReplaceAll("a", "");
    drawAntarctica = true;
  }

  if(drawAntarctica){
    fAntarctica = getAntarctica();
    fAntarctica->Draw("colz");

    // now draw on top of Antarctica
    opt += "same";
  }

  // now call the regular TGraph::Draw()
  TGraph::Draw(opt);

  // now for the hack...
  // does this redirect scaling etc.?
  fHistogram = (TH1F*) fAntarctica;

}



void TGraphAntarctica::init(){

  const int defaultCoarsenessFactor = 50;

  cf = defaultCoarsenessFactor;
  lastCf = cf;

  dataSet = RampdemReader::surface;
  lastDataSet = dataSet;

  fAntarctica = NULL;

  doneConversion = false;
  convertArrays();
}



TProfile2D* TGraphAntarctica::getAntarctica(){

  if(!fAntarctica){
    fAntarctica = RampdemReader::getMap(dataSet, cf);
  }
  else if(cf!=lastCf || dataSet!=lastDataSet){
    delete fAntarctica;
    fAntarctica = RampdemReader::getMap(dataSet, cf);
  }

  lastDataSet = dataSet;
  lastCf = cf;

  fAntarctica->SetDirectory(0);

  return fAntarctica;
}


void TGraphAntarctica::convertArrays(){
  if(!doneConversion){
    for(Int_t i=0; i < GetN(); i++){
      RampdemReader::LonLatToEastingNorthing(fX[i], fY[i], fX[i], fY[i]);
    }
  }
  doneConversion = true;
}
