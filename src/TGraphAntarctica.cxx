#include "TGraphAntarctica.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "Adu5Pat.h"
#include "AnitaDataset.h"
#include "TVirtualPad.h"

ClassImp(TGraphAntarctica)


/** 
 * Replace all characters that make it hard to interact with on the root prompt
 * 
 * @param name the name to sanitize
 * 
 * @return The sanitized name
 */
TString makeSanitizedName(const char* name){
  TString sanitizedName(name);

  // try to remove anything that would confuse getting this name
  // from a file by typing the name on the command line
  sanitizedName.ReplaceAll(" ", "_");
  sanitizedName.ReplaceAll("-", "_");
  sanitizedName.ReplaceAll("(", "");
  sanitizedName.ReplaceAll(")", "");
  sanitizedName.ReplaceAll(")", "");
  sanitizedName.ReplaceAll(";", "");
  sanitizedName.ReplaceAll(".", "");
  return sanitizedName;
}



/** 
 * Set the ith point of the Graph to a longitude/latitude
 * Handles the conversion to easting/northing internally
 * 
 * @param i is the graph point
 * @param lon is the longitude
 * @param lat is the latitude
 */
void TGraphAntarctica::SetPoint(Int_t i, Double_t lon, Double_t lat){
  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  TGraph::SetPoint(i, easting, northing);
}



/** 
 * Set the ith point to the position held in the AntarcticCoord
 * 
 * @param i is the point to set 
 * @param coord is the position (is internally converted to WGS84)
 */
void TGraphAntarctica::SetPoint(Int_t i, const AntarcticCoord& coord){
  AntarcticCoord stereo = coord.as(AntarcticCoord::STEREOGRAPHIC);
  TGraph::SetPoint(i, stereo.x, stereo.y);
}




/** 
 * Construct a TGraph antarctica from run firstRun to lastRun (inclusive)
 * This is mostly for plotting purposes, otherwise 
 * 
 * @param firstRun is the first run
 * @param lastRun is the last run
 * @param pointEvery 
 * @param quiet default true
 *
 * @return the newly constructed TGraphAntarctica
 */
TGraphAntarctica* TGraphAntarctica::makeGpsGraph(int firstRun, int lastRun, int gpsTreeStride,bool quiet){
  if (!quiet) std::cout << "makeGpsGraph(): starting..." << std::endl;

  // handle default tree stride
  gpsTreeStride = gpsTreeStride <= 0 ? defaultGpsTreeStride : gpsTreeStride;

  TGraphAntarctica* gr = new TGraphAntarctica();

  for(int run=firstRun; run<=lastRun; run++){
    if (AnitaVersion::get() == 3) {
      if (run > 256 && run < 264) {
	std::cout << "makeGpsGraph(): In ANITA3 runs 257 through 263 are broken, skipping to 264..." << std::endl;
	run = 264;
      }
    }

    AnitaDataset d(run);
    if (!quiet) std::cout << "makeGpsGraph(): starting run" << run << " - d.N()=" << d.N() << std::endl;
    for(int entry=0; entry < d.N(); entry+=gpsTreeStride){
      if (!quiet) std::cout << "makeGpsGraph(): run:" << run << " entry:" << entry << std::endl;
      d.getEntry(entry);
      Adu5Pat* pat = d.gps();

      gr->SetPoint(gr->GetN(), pat->longitude, pat->latitude);

    }
  }
  return gr;
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



/** 
 * Constructor for BaseList::path
 * 
 * @param p is the const reference to the path
 * @param interpSeconds is a time to interpolate between point, if 0 the raw data is used and no interpolation is done (default=0)
 */
TGraphAntarctica::TGraphAntarctica(const BaseList::path& p, UInt_t interpSeconds){
  init();
  if(interpSeconds == 0){
    for(UInt_t i=0; i < p.ps.size(); i++){
      SetPoint(i, p.ps.at(i));
    }
  }
  else{
    UInt_t time = p.ts.at(0);
    while(time < p.ts.back()){
      SetPoint(GetN(), p.getPosition(time));
      time += interpSeconds;      
    }
  }
  SetEditable(false);
  SetNameTitle(makeSanitizedName(p.getName()), p.getName());
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
  fBackground = NULL;
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



TGraphAntarctica::TGraphAntarctica(const BaseList::base& b){
  init();
  SetPoint(0, b.getPosition(0));
  SetNameTitle(makeSanitizedName(b.getName()), b.getName());
  SetEditable(false);  
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



TAxis* TGraphAntarctica::GetXaxis(){
  AntarcticaBackground* b = getBackground();
  if(b){
    return b->GetXaxis();
  }
  else{
    return fHistogram->GetXaxis();
  }
}

TAxis* TGraphAntarctica::GetYaxis(){
  AntarcticaBackground* b = getBackground();
  if(b){
    return b->GetYaxis();
  }
  else{
    return fHistogram->GetYaxis();
  }
}
