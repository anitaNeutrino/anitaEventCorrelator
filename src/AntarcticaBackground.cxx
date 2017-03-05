#include "AntarcticaBackground.h"
#include "TGraphAntarctica.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"


ClassImp(AntarcticaBackground)




static int numAntarcticaBackgrounds = 0;


/**
 * Helper function to generate the an instance of the Antarctica Background.
 * This is what you should use to make an instance in basically all circumstances.
 *
 * @param dataSet
 * @param coarseness
 *
 * @return
 */
AntarcticaBackground* AntarcticaBackground::generate(RampdemReader::dataSet dataSet, Int_t coarseness){

  Int_t nx, ny;
  RampdemReader::getNumXY(nx, ny, dataSet);
  Double_t xMin, xMax, yMin, yMax;
  RampdemReader::getMapCoordinates(xMin, yMin, xMax, yMax, dataSet);

  nx /= coarseness;
  ny /= coarseness;

  TString name = TString::Format("fAntarctica%d", numAntarcticaBackgrounds);
  AntarcticaBackground* ab = new AntarcticaBackground(dataSet, coarseness, name, "", nx, xMin, xMax, ny, yMin, yMax, "");
  numAntarcticaBackgrounds++;
  return ab;
}


AntarcticaBackground::AntarcticaBackground(){
  std::cerr << "Nope" << std::endl;
}





/**
 * The useful constructor called by AntarcticaBackground::generate
 *
 * @param dataSet
 * @param coarseness
 * @param name
 * @param title
 * @param nx
 * @param xlow
 * @param xup
 * @param ny
 * @param ylow
 * @param yup
 * @param option
 */
AntarcticaBackground::AntarcticaBackground(RampdemReader::dataSet dataSet, Int_t coarseness,
					   const char *name,const char *title,
					   Int_t nx,Double_t xlow,Double_t xup,
					   Int_t ny,Double_t ylow,Double_t yup,
					   Option_t *option) :
  TProfile2D(name, title, nx, xlow, xup, ny, ylow, yup, option){

  fDataSet = dataSet;
  fCoarseness = coarseness;
  std::cerr << "Yep" << std::endl;

  // force update to read in data...
  lastDataSet = fDataSet;
  lastCoarseness = -1;

  update();
}



/**
 * Default Constructor
 */
// AntarcticaBackground::AntarcticaBackground(RampdemReader::dataSet ds, Int_t coarseness){

//   fName = "TestyMcTesterson";

//   fCoarseness = coarseness;
//   lastCoarseness = fCoarseness;

//   fDrawLonLatGrids = true;
//   fLonLatGridPoints = 36000;
//   lastLonLatGridPoints = fLonLatGridPoints;

//   fDataSet = ds;
//   lastDataSet = fDataSet;

//   RampdemReader::getNumXY(fNumX, fNumY);

//   update();
// }


/**
 * Default Destructor
 */
// AntarcticaBackground::~AntarcticaBackground(){
//   // if(fAntarctica){
//   //   delete fAntarctica;
//   // }
// }



/**
 * Update the map of Antarctica as different options are selected
 */
void AntarcticaBackground::update(){


  Int_t nx, ny;
  RampdemReader::getNumXY(nx, ny, fDataSet);
  Double_t xMin, xMax, yMin, yMax;
  RampdemReader::getMapCoordinates(xMin, yMin, xMax, yMax, fDataSet);

  nx /= fCoarseness;
  ny /= fCoarseness;


  if(fCoarseness!=lastCoarseness || fDataSet!=lastDataSet){

    // emptry the TProfile...
    fBinEntries.Reset();
    for(int by=0; by <= GetNbinsY() + 1; by++){
      for(int bx=0; bx <= GetNbinsX() + 1; bx++){
	SetBinContent(bx, by, 0);
      }
    }

    // change the dimensions?
    SetBins(nx, xMin, xMax, ny, yMin, yMax);

    // insert new data
    RampdemReader::fillThisHist(this, fDataSet);
  }

  // book keeping and prettification
  // fAntarctica->SetName("fAntarctica");
  GetXaxis()->SetNdivisions(0, kFALSE);
  GetYaxis()->SetNdivisions(0, kFALSE);

  lastDataSet = fDataSet;
  lastCoarseness = fCoarseness;

  SetDirectory(0);

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
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", coarsenesss must be >= 1." << std::endl;
    coarseness = 1;
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



// void AntarcticaBackground::setColAxisTitle(){

//   TObjArray* tokens = fTitle.Tokenize(";");

//   TString newTitle = tokens->GetEntries() > 0 ? ((TObjString*) tokens->At(0))->GetString() : "";

//   // pad for x and y axes (not drawn)
//   newTitle += "; ; ;";

//   newTitle += TString::Format("%s", RampdemReader::dataSetToAxisTitle(fDataSet));

// }


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


// void AntarcticaBackground::makePrettyPalette(){
//   gPad->Modified();
//   gPad->Update();
//   TPaletteAxis *palette = (TPaletteAxis*) GetListOfFunctions()->FindObject("palette");
//   palette->SetX1NDC(0.03);
//   palette->SetX2NDC(0.06);
//   palette->SetY1NDC(0.03);
//   palette->SetY2NDC(0.16);
//   palette->SetTitleSize(0.001);
//   palette->SetTitleOffset(0.1);
//   gPad->Modified();
//   gPad->Update();

// }


// void AntarcticaBackground::setPadMargins(){
//   gPad->SetTopMargin(0.05);
//   gPad->SetBottomMargin(0.05);
//   gPad->SetLeftMargin(0.05);
//   gPad->SetRightMargin(0.05);
//   gPad->SetFrameLineColor(0);
//   gPad->SetFrameLineWidth(0);
//   gPad->SetFrameBorderSize(0);
// }















void AntarcticaBackground::SetRampdemDataSet(bool useRampdemDataSet){
  SetDataSet(RampdemReader::rampdem);
}

Bool_t AntarcticaBackground::GetRampdemDataSet(){
  return fDataSet == RampdemReader::rampdem;
}

void AntarcticaBackground::SetBedDataSet(bool useBedDataSet){
  SetDataSet(RampdemReader::bed);
};

Bool_t AntarcticaBackground::GetBedDataSet(){
  return fDataSet == RampdemReader::bed;
}

void AntarcticaBackground::SetIcemaskDataSet(bool useIcemaskDataSet){
  SetDataSet(RampdemReader::icemask_grounded_and_shelves);
}

Bool_t AntarcticaBackground::GetIcemaskDataSet(){
  return fDataSet == RampdemReader::icemask_grounded_and_shelves;
}

void AntarcticaBackground::SetSurfaceDataSet(bool useSurfaceDataSet){
  SetDataSet(RampdemReader::surface);
}

Bool_t AntarcticaBackground::GetSurfaceDataSet(){
  return fDataSet == RampdemReader::surface;
}

void AntarcticaBackground::SetThicknessDataSet(bool useThicknessDataSet){
  SetDataSet(RampdemReader::thickness);
}

Bool_t AntarcticaBackground::GetThicknessDataSet(){
  return fDataSet == RampdemReader::thickness;
}
