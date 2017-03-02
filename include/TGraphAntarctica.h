#include "TGraph.h"
#include "RampdemReader.h"
#include "TVirtualGraphPainter.h"






  // here I hope to overload all the functions which can affect the X, Y arrays inside a TGraph
  // and add a step to convert X=lon, Y=lat, to X=easting, Y=northing for pretty plotting
  // I may even overload the draw function to add an option to draw an Antarctic Histogram in the background

class TGraphAntarctica : public TGraph {

private:

  Int_t cf; // map coarseness factor
  RampdemReader::dataSet dataSet;
  Bool_t doneConversion;
  void convertArrays(){
    if(!doneConversion){
      for(Int_t i=0; i < GetN(); i++){
	RampdemReader::LonLatToEastingNorthing(fX[i], fY[i], fX[i], fY[i]);
      }
    }
    doneConversion = true;
  }

public:

  TGraphAntarctica(Int_t n, const Int_t *x, const Int_t *y) : TGraph(n, x, y) {
    convertArrays();
    dataSet = RampdemReader::surface;
  }
  TGraphAntarctica(Int_t n, const Float_t *x, const Float_t *y) : TGraph(n, x, y) {
    convertArrays();
    dataSet = RampdemReader::surface;
  };
  TGraphAntarctica(Int_t n, const Double_t *x, const Double_t *y) : TGraph(n, x, y) {
    convertArrays();
    dataSet = RampdemReader::surface;
  }

  explicit TGraphAntarctica(const TGraph &gr) : TGraph(gr) {
    convertArrays();
    dataSet = RampdemReader::surface;
  }


  // TGraphAntarctica& operator=(const TGraph& gr){
  //   TGraph::operator=(gr);
  //   convertArrays();
  // }

  // TGraphAntarctica& operator=(const TGraphAntarctica& gr){
  //   TGraph::operator=(gr);
  // }

  explicit TGraphAntarctica(const TVectorF &vx, const TVectorF &vy) : TGraph(vx, vy) {
    convertArrays();
    dataSet = RampdemReader::surface;
  }

  explicit TGraphAntarctica(const TVectorD &vx, const TVectorD &vy) : TGraph(vx, vy) {
    convertArrays();
    dataSet = RampdemReader::surface;
  }
  explicit TGraphAntarctica(const TH1 *h) : TGraph(h){
    convertArrays();
    dataSet = RampdemReader::surface;
  }

  explicit TGraphAntarctica(const TF1 *f, Option_t *option="") : TGraph(f, option){
    convertArrays();
    dataSet = RampdemReader::surface;
  }

  explicit TGraphAntarctica(const char *filename, const char *format="%lg %lg", Option_t *option="") : TGraph(filename, format, option) {
    convertArrays();
  }

  void SetPoint(Int_t i, Double_t lon, Double_t lat){
    Double_t easting, northing;
    RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
    TGraph::SetPoint(i, easting, northing);
  }

  Double_t* GetEasting(){
    return GetX();
  }
  Double_t* GetNorthing(){
    return GetY();
  }

  TH1F *GetHistogram() const; // where the magic happens

};
