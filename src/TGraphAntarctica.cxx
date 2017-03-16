#include "TGraphAntarctica.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "Adu5Pat.h"

ClassImp(TGraphAntarctica)


// Don't want this to be part of the class when we start writing things to files.
static AntarcticaBackground* fBackground = NULL;


void TGraphAntarctica::SetPoint(Int_t i, Double_t lon, Double_t lat){
  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  TGraph::SetPoint(i, easting, northing);
}








TGraphAntarctica::TGraphAntarctica(TChain* chain, TString lonSelector, TString latSelector, TCut cut) : TGraph(){

  TString command = lonSelector + ":" + latSelector;

  const int n = chain->Draw(command, cut, "goff"); // "goff" means graphics off
  Set(n);
  Double_t* lons = chain->GetV1();
  Double_t* lats = chain->GetV2();

  for(int i=0; i < n; i++){
    fX[i] = lons[i];
    fY[i] = lats[i];
    // std::cout << i << "\t" << fX[i] << "\t" << fY[i] << std::endl;
  }
  init();
}

TGraphAntarctica::TGraphAntarctica(TTree* tree, TString lonSelector, TString latSelector, TCut cut) : TGraph(){

  TString command = lonSelector + ":" + latSelector;

  const int n = tree->Draw(command, cut, "goff"); // "goff" means graphics off
  Set(n);
  Double_t* lons = tree->GetV1();
  Double_t* lats = tree->GetV2();

  for(int i=0; i < n; i++){
    fX[i] = lons[i];
    fY[i] = lats[i];
    // std::cout << i << "\t" << fX[i] << "\t" << fY[i] << std::endl;
  }
  init();
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
    // std::cout << "prims " << prims->GetEntries() << std::endl;
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

    const char* extraText = fTitle.Length() > 0 ? fTitle.Data() : NULL;
    fBackground->updateToolTip(event, x, y, extraText);
  }

  TGraph::ExecuteEvent(event, x, y);
}







// enum EColorPalette {kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
//                     kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
//                     kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
//                     kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
//                     kAlpine=63,           kAquamarine=64,   kArmy=65,
//                     kAtlantic=66,         kAurora=67,       kAvocado=68,
//                     kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
//                     kBrownCyan=72,        kCMYK=73,         kCandy=74,
//                     kCherry=75,           kCoffee=76,       kDarkRainBow=77,
//                     kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
//                     kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
//                     kGreenPink=84,        kIsland=85,       kLake=86,
//                     kLightTemperature=87, kLightTerrain=88, kMint=89,
//                     kNeon=90,             kPastel=91,       kPearl=92,
//                     kPigeon=93,           kPlum=94,         kRedBlue=95,
//                     kRose=96,             kRust=97,         kSandyTerrain=98,
//                     kSienna=99,           kSolar=100,       kSouthWest=101,
//                     kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
//                     kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
//                     kWaterMelon=108,      kCool=109,        kCopper=110,
//                     kGistEarth=111,       kViridis=112};
