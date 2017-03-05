#include "AntarcticaBackground.h"
#include "TGraphAntarctica.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"


ClassImp(AntarcticaBackground)



// simple global variable to increment integer in name so we can have as manny of these background as we like.
static int numAntarcticaBackgrounds = 0;


AntarcticaBackground::AntarcticaBackground(RampdemReader::dataSet dataSet, Int_t coarseness)
  : TProfile2D() {
  init(dataSet, coarseness);
}


void AntarcticaBackground::init(RampdemReader::dataSet dataSet, Int_t coarseness){
  SetDirectory(0);
  fName = TString::Format("fAntarctica%d", numAntarcticaBackgrounds);
  numAntarcticaBackgrounds++;
  fDataSet = dataSet;
  fCoarseness = coarseness;

  lastDataSet = fDataSet;
  lastCoarseness = -1;  // force update() to read in data by setting impossible coarseness
  update();

}



/**
 * Update the map of Antarctica as different options are selected
 */
void AntarcticaBackground::update(){


  if(fCoarseness!=lastCoarseness || fDataSet!=lastDataSet){

    // Here I save the current drawn axis range so I can set it to be the same after updating contents
    Int_t firstX = fXaxis.GetFirst();
    Int_t lastX = fXaxis.GetLast();
    Int_t firstY = fYaxis.GetFirst();
    Int_t lastY = fYaxis.GetLast();

    Double_t lowX = fXaxis.GetBinLowEdge(firstX);
    Double_t highX = fXaxis.GetBinUpEdge(lastX);
    Double_t lowY = fXaxis.GetBinLowEdge(firstY);
    Double_t highY = fXaxis.GetBinUpEdge(lastY);

    std::cout << firstX << "\t" << lastX << "\t" << lowX << "\t" << highX << std::endl;


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

    if(gPad){

      // now set the viewing range to the same as before if we have a gPad instance
      fXaxis.SetRangeUser(lowX, highX);
      fYaxis.SetRangeUser(lowY, highY);

      // and update the z-axis title if needed
      prettifyPalette();
    }

  }

  // prettification
  GetXaxis()->SetNdivisions(0, kFALSE);
  GetYaxis()->SetNdivisions(0, kFALSE);

  lastDataSet = fDataSet;
  lastCoarseness = fCoarseness;

  // if(alreadyDrawn){

  //   setPadMargins();

  // if(gPad){
  //   TList* prims = gPad->GetListOfPrimitives();
  //   //   fAntarctica->SetOption("col2azfbbb");

  //   //   prims->AddBefore(this, fAntarctica);
  //   //   // fHistogram = (TH1F*) fAntarctica;

  //   if(fDrawLonLatGrids){
  //     for(UInt_t grInd=0; grInd < grLonLatGrids.size(); grInd++){
  // 	TGraph* gr = grLonLatGrids.at(grInd);
  // 	// gr->SetOption("lsame");
  // 	prims->AddBefore(this, gr);
  //     }
  //   }
  //   else{
  //     for(UInt_t grInd=0; grInd < grLonLatGrids.size(); grInd++){
  // 	TGraph* gr = grLonLatGrids.at(grInd);
  // 	// gr->SetOption("lsame");
  // 	prims->RecursiveRemove(gr);
  //     }
  //   }
  // }
  // }
  // fHistogram = (TH1F*) fAntarctica;

  // return fAntarctica;

}



// /**
//  * @brief Draws the histogram
//  *
//  * @param option is the draw option, which is passed on to the histogram
//  */
// void AntarcticaBackground::Draw(Option_t* option){

//   TString opt = option;
//   opt.ToLower();

//   // If there is no pad or an empty pad the "same" option is ignored.
//   if (!gPad) {
//     gROOT->MakeDefCanvas();
//   }

//   setPadMargins();

//   AppendPad(opt.Data());


// }


// void AntarcticaBackground::Paint(Option_t* opt){

//   fAntarctica->Paint(opt);
//   makePrettyPalette();
//   setPadMargins();
// }






Int_t AntarcticaBackground::GetCoarseness(){
  return fCoarseness;
}


void AntarcticaBackground::SetCoarseness(Int_t coarseness){
  if(coarseness < 1){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", coarsenesss must be >= 1. Setting coareness = "
	      << defaultCoarseness << std::endl;
    coarseness = defaultCoarseness;
  }
  fCoarseness = coarseness;
  update();
}


RampdemReader::dataSet AntarcticaBackground::GetDataSet(){
  return fDataSet;
}



void AntarcticaBackground::SetDataSet(RampdemReader::dataSet dataSet){

  // this is all that's really needed...
  fDataSet = dataSet;
  update();
}





// void AntarcticaBackground::SetDrawLonLatGrids(Bool_t drawLonLatGrids){
//   fDrawLonLatGrids = drawLonLatGrids;
//   update();
// }

// Bool_t AntarcticaBackground::GetDrawLonLatGrids(){
//   return fDrawLonLatGrids;
// }



// void AntarcticaBackground::makeLonLatGrids(){

//   const int minLat = -85;
//   const int maxLat = -60;
//   const int deltaLon = 45; // degrees
//   const int deltaLat = 5; // degrees

//   if(fLonLatGridPoints != lastLonLatGridPoints){
//     deleteLonLatGrids();

//     // make circles of constant latitude
//     for(int lat = minLat; lat<= maxLat; lat += deltaLat){
//       TGraph* gr = new TGraph();
//       gr->SetLineColor(kGray);
//       const Double_t deltaLon = 360./fLonLatGridPoints;
//       for(int i=0; i < fLonLatGridPoints; i++){
// 	Double_t theLat = lat;
// 	Double_t theLon = i*deltaLon;
// 	Double_t easting, northing;
// 	RampdemReader::LonLatToEastingNorthing(theLon, theLat, easting, northing);
// 	gr->SetPoint(gr->GetN(), easting, northing);
// 	// std::cout << gr << "\t" << gr->GetN() << "\t" << easting << "\t" << northing << std::endl;
//       }
//       grLonLatGrids.push_back(gr);
//     }

//     // make lines of constant longitude
//     for(int lon = 0; lon < 360; lon+= deltaLon){
//       TGraph* gr = new TGraph();
//       gr->SetLineColor(kGray);
//       const Double_t deltaLat = double(maxLat - -90)/fLonLatGridPoints;
//       for(int i=0; i < fLonLatGridPoints; i++){
// 	Double_t theLat = -90 + deltaLat*i;
// 	Double_t theLon = lon;
// 	Double_t easting, northing;
// 	RampdemReader::LonLatToEastingNorthing(theLon, theLat, easting, northing);
// 	// std::cout << gr << "\t" << gr->GetN() << "\t" << easting << "\t" << northing << "\t" << theLat << "\t" << theLon << std::endl;
// 	gr->SetPoint(gr->GetN(), easting, northing);
//       }
//       grLonLatGrids.push_back(gr);
//     }

//     lastLonLatGridPoints = fLonLatGridPoints;
//   }

// }




// void AntarcticaBackground::deleteLonLatGrids(){
//   while(grLonLatGrids.size() > 0){
//     TGraph* gr = grLonLatGrids.back();
//     delete gr;
//     grLonLatGrids.pop_back();
//   }

// }










/**
 * Draw function for Antarctica Background which also prettifies the pad/canvas.
 *
 * @param opt is the draw option (default is colz)
 */
void AntarcticaBackground::Draw(Option_t* opt){
  TProfile2D::Draw(opt);
  setPadMargins();
  prettifyPalette();
  fXaxis.SetAxisColor(kWhite);
  fYaxis.SetAxisColor(kWhite);
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
    std::cout << zAxis->GetTitleOffset() << std::endl;
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
