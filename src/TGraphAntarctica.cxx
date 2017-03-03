#include "TGraphAntarctica.h"
#include "TVirtualPad.h"

ClassImp(TGraphAntarctica)


TGraphAntarctica::~TGraphAntarctica(){
  if(fAntarctica){
    delete fAntarctica;
  }
  deleteLonLatGrids();
}


void TGraphAntarctica::SetDrawLonLatGrids(Bool_t drawLonLatGrids){
  fDrawLonLatGrids = drawLonLatGrids;
  fAntarctica = getAntarctica();
}

Bool_t TGraphAntarctica::GetDrawLonLatGrids(){
  return fDrawLonLatGrids;
}



void TGraphAntarctica::SetPoint(Int_t i, Double_t lon, Double_t lat){
  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  TGraph::SetPoint(i, easting, northing);
}


Int_t TGraphAntarctica::GetCoarseness(){
  return fCoarseness;
}

void TGraphAntarctica::SetCoarseness(Int_t coarseness){
  if(coarseness < 1){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", coarsenesss must be >= 1." << std::endl;
    coarseness = 1;
  }
  fCoarseness = coarseness;
  fAntarctica = getAntarctica();
}


RampdemReader::dataSet TGraphAntarctica::GetDataSet(){
  return fDataSet;
}



void TGraphAntarctica::SetDataSet(RampdemReader::dataSet dataSet){

  // this is all that's really needed...
  fDataSet = dataSet;
  fAntarctica = getAntarctica();
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

  if(fDrawLonLatGrids){
    makeLonLatGrids();
    for(UInt_t grInd=0; grInd < grLonLatGrids.size(); grInd++){
      grLonLatGrids.at(grInd)->Draw("lsame");
    }
  }

  // now call the regular TGraph::Draw()
  TGraph::Draw(opt);

  // now for the hack...
  // does this redirect scaling etc.?
  fHistogram = (TH1F*) fAntarctica;

  alreadyDrawn = true;

}



void TGraphAntarctica::init(){

  const int defaultCoarsenessFactor = 10;

  fCoarseness = defaultCoarsenessFactor;
  lastCoarseness = fCoarseness;

  fAntarctica = NULL;
  alreadyDrawn = false;

  fDrawLonLatGrids = true;
  fLonLatGridPoints = 36000;
  lastLonLatGridPoints = fLonLatGridPoints;

  doneConversion = false;
  convertArrays();


  // default data set is set in this map...
  // make sure only one of these is true!
  dataSetBools[RampdemReader::rampdem] = false;
  dataSetBools[RampdemReader::bed] = false;
  dataSetBools[RampdemReader::icemask_grounded_and_shelves] = false;
  dataSetBools[RampdemReader::surface] = true;
  dataSetBools[RampdemReader::thickness] = false;

  for(std::map<RampdemReader::dataSet, Bool_t>::iterator it=dataSetBools.begin(); it!=dataSetBools.end(); ++it){
    if(it->second == true){
      fDataSet = it->first;
      lastDataSet = fDataSet;
      break;
    }
  }


}




TProfile2D* TGraphAntarctica::getAntarctica(){

  if(!fAntarctica){
    fAntarctica = RampdemReader::getMap(fDataSet, fCoarseness);
  }
  else if(fCoarseness!=lastCoarseness || fDataSet!=lastDataSet){
    delete fAntarctica;
    fAntarctica = RampdemReader::getMap(fDataSet, fCoarseness);
  }

  lastDataSet = fDataSet;
  lastCoarseness = fCoarseness;

  fAntarctica->SetDirectory(0);

  if(alreadyDrawn){
    TList* prims = gPad->GetListOfPrimitives();
    fAntarctica->SetOption("colz");
    prims->AddBefore(this, fAntarctica);
    fHistogram = (TH1F*) fAntarctica;

    if(fDrawLonLatGrids){
      for(UInt_t grInd=0; grInd < grLonLatGrids.size(); grInd++){
	TGraph* gr = grLonLatGrids.at(grInd);
	// gr->SetOption("lsame");
	prims->AddBefore(this, gr);
      }
    }
    else{
      for(UInt_t grInd=0; grInd < grLonLatGrids.size(); grInd++){
	TGraph* gr = grLonLatGrids.at(grInd);
	// gr->SetOption("lsame");
	prims->RecursiveRemove(gr);
      }
    }
  }
  fHistogram = (TH1F*) fAntarctica;

  return fAntarctica;
}


void TGraphAntarctica::makeLonLatGrids(){

  const int minLat = -85;
  const int maxLat = -60;
  const int deltaLon = 45; // degrees
  const int deltaLat = 5; // degrees

  if(!alreadyDrawn || fLonLatGridPoints != lastLonLatGridPoints){
    deleteLonLatGrids();

    // make circles of constant latitude
    for(int lat = minLat; lat<= maxLat; lat += deltaLat){
      TGraph* gr = new TGraph();
      gr->SetLineColor(kGray);
      const Double_t deltaLon = 360./fLonLatGridPoints;
      for(int i=0; i < fLonLatGridPoints; i++){
	Double_t theLat = lat;
	Double_t theLon = i*deltaLon;
	Double_t easting, northing;
	RampdemReader::LonLatToEastingNorthing(theLon, theLat, easting, northing);
	gr->SetPoint(gr->GetN(), easting, northing);
	// std::cout << gr << "\t" << gr->GetN() << "\t" << easting << "\t" << northing << std::endl;
      }
      grLonLatGrids.push_back(gr);
    }

    // make lines of constant longitude
    for(int lon = 0; lon < 360; lon+= deltaLon){
      TGraph* gr = new TGraph();
      gr->SetLineColor(kGray);
      const Double_t deltaLat = double(maxLat - -90)/fLonLatGridPoints;
      for(int i=0; i < fLonLatGridPoints; i++){
	Double_t theLat = -90 + deltaLat*i;
	Double_t theLon = lon;
	Double_t easting, northing;
	RampdemReader::LonLatToEastingNorthing(theLon, theLat, easting, northing);
	// std::cout << gr << "\t" << gr->GetN() << "\t" << easting << "\t" << northing << "\t" << theLat << "\t" << theLon << std::endl;
	gr->SetPoint(gr->GetN(), easting, northing);
      }
      grLonLatGrids.push_back(gr);
    }

    lastLonLatGridPoints = fLonLatGridPoints;
  }

}


void TGraphAntarctica::convertArrays(){
  if(!doneConversion){
    for(Int_t i=0; i < GetN(); i++){
      RampdemReader::LonLatToEastingNorthing(fX[i], fY[i], fX[i], fY[i]);
    }
  }

  doneConversion = true;
}


void TGraphAntarctica::deleteLonLatGrids(){
  while(grLonLatGrids.size() > 0){
    TGraph* gr = grLonLatGrids.back();
    delete gr;
    grLonLatGrids.pop_back();
  }

}






void TGraphAntarctica::SetRampdemDataSet(bool useRampdemDataSet){
  SetDataSet(RampdemReader::rampdem);
}

Bool_t TGraphAntarctica::GetRampdemDataSet(){
  return fDataSet == RampdemReader::rampdem;
}

void TGraphAntarctica::SetBedDataSet(bool useBedDataSet){
  SetDataSet(RampdemReader::bed);
};

Bool_t TGraphAntarctica::GetBedDataSet(){
  return fDataSet == RampdemReader::bed;
}

void TGraphAntarctica::SetIcemaskDataSet(bool useIcemaskDataSet){
  SetDataSet(RampdemReader::icemask_grounded_and_shelves);
}

Bool_t TGraphAntarctica::GetIcemaskDataSet(){
  return fDataSet == RampdemReader::icemask_grounded_and_shelves;
}

void TGraphAntarctica::SetSurfaceDataSet(bool useSurfaceDataSet){
  SetDataSet(RampdemReader::surface);
}

Bool_t TGraphAntarctica::GetSurfaceDataSet(){
  return fDataSet == RampdemReader::surface;
}

void TGraphAntarctica::SetThicknessDataSet(bool useThicknessDataSet){
  SetDataSet(RampdemReader::thickness);
}

Bool_t TGraphAntarctica::GetThicknessDataSet(){
  return fDataSet == RampdemReader::thickness;
}
