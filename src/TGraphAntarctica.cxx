#include "TGraphAntarctica.h"
#include "TVirtualPad.h"


void TGraphAntarctica::SetPoint(Int_t i, Double_t lon, Double_t lat){
  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  TGraph::SetPoint(i, easting, northing);
}


Int_t TGraphAntarctica::GetCoarseness(){
  return cf;
}

void TGraphAntarctica::SetCoarseness(Int_t newCf){
  if(newCf < 1){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", coarsenesss must be >= 1." << std::endl;
    newCf = 1;
  }
  cf = newCf;
}


RampdemReader::dataSet TGraphAntarctica::GetDataSet(){
  return dataSet;
}


void TGraphAntarctica::SetDataSet(RampdemReader::dataSet newDataSet){
  dataSet = newDataSet;
}


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

  const int defaultCoarsenessFactor = 10;

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
