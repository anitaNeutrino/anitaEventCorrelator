#include "TGraphAntarcticaa.h"
#include "TVirtualPad.h"



TH1F* TGraphAntarctica::GetHistogram() const{

  // Returns a pointer to the histogram used to draw the axis
  // Takes into account the two following cases.
  //    1- option 'A' was specified in TGraph::Draw. Return fHistogram
  //    2- user had called TPad::DrawFrame. return pointer to hframe histogram

  Double_t rwxmin, rwxmax, rwymin, rwymax, maximum, minimum, dx, dy;
  Double_t uxmin, uxmax;

  ComputeRange(rwxmin, rwymin, rwxmax, rwymax);  //this is redefined in TGraphErrors

  // (if fHistogram exist) && (if the log scale is on) &&
  // (if the computed range minimum is > 0) && (if the fHistogram minimum is zero)
  // then it means fHistogram limits have been computed in linear scale
  // therefore they might be too strict and cut some points. In that case the
  // fHistogram limits should be recomputed ie: the existing fHistogram
  // should not be returned.
  TH1F *historg = 0;
  if (fHistogram) {
    if (gPad && gPad->GetLogx()) {
      if (rwxmin <= 0 || fHistogram->GetXaxis()->GetXmin() != 0) return fHistogram;
    } else if (gPad && gPad->GetLogy()) {
      if (rwymin <= 0 || fHistogram->GetMinimum() != 0) return fHistogram;
    } else {
      return fHistogram;
    }
    historg = fHistogram;
  }

  if (rwxmin == rwxmax) rwxmax += 1.;
  if (rwymin == rwymax) rwymax += 1.;
  dx = 0.1 * (rwxmax - rwxmin);
  dy = 0.1 * (rwymax - rwymin);
  uxmin    = rwxmin - dx;
  uxmax    = rwxmax + dx;
  minimum  = rwymin - dy;
  maximum  = rwymax + dy;
  if (fMinimum != -1111) minimum = fMinimum;
  if (fMaximum != -1111) maximum = fMaximum;

  // the graph is created with at least as many channels as there are points
  // to permit zooming on the full range
  if (uxmin < 0 && rwxmin >= 0) {
    if (gPad && gPad->GetLogx()) uxmin = 0.9 * rwxmin;
    else                 uxmin = 0;
  }
  if (uxmax > 0 && rwxmax <= 0) {
    if (gPad && gPad->GetLogx()) uxmax = 1.1 * rwxmax;
    else                 uxmax = 0;
  }
  if (minimum < 0 && rwymin >= 0) {
    if (gPad && gPad->GetLogy()) minimum = 0.9 * rwymin;
    else                minimum = 0;
  }
  if (minimum <= 0 && gPad && gPad->GetLogy()) minimum = 0.001 * maximum;
  if (uxmin <= 0 && gPad && gPad->GetLogx()) {
    if (uxmax > 1000) uxmin = 1;
    else              uxmin = 0.001 * uxmax;
  }

  rwxmin = uxmin;
  rwxmax = uxmax;
  Int_t npt = 100;
  if (fNpoints > npt) npt = fNpoints;
  const char *gname = GetName();
  if (!gname[0]) gname = "GraphAntarctica";
  // ((TGraph*)this)->fHistogram = new TH1F(gname, GetTitle(), npt, rwxmin, rwxmax);
  ((TGraphAntarctica*)this)->fHistogram = RampdemReader::getMap(dataSet, cf);

  if (!fHistogram) return 0;
  fHistogram->SetMinimum(minimum);
  fHistogram->SetBit(TH1::kNoStats);
  fHistogram->SetMaximum(maximum);
  fHistogram->GetYaxis()->SetLimits(minimum, maximum);
  fHistogram->SetDirectory(0);
  // Restore the axis attributes if needed
  if (historg) {
    fHistogram->GetXaxis()->SetTitle(historg->GetXaxis()->GetTitle());
    fHistogram->GetXaxis()->CenterTitle(historg->GetXaxis()->GetCenterTitle());
    fHistogram->GetXaxis()->RotateTitle(historg->GetXaxis()->GetRotateTitle());
    fHistogram->GetXaxis()->SetNoExponent(historg->GetXaxis()->GetNoExponent());
    fHistogram->GetXaxis()->SetNdivisions(historg->GetXaxis()->GetNdivisions());
    fHistogram->GetXaxis()->SetLabelFont(historg->GetXaxis()->GetLabelFont());
    fHistogram->GetXaxis()->SetLabelOffset(historg->GetXaxis()->GetLabelOffset());
    fHistogram->GetXaxis()->SetLabelSize(historg->GetXaxis()->GetLabelSize());
    fHistogram->GetXaxis()->SetTitleSize(historg->GetXaxis()->GetTitleSize());
    fHistogram->GetXaxis()->SetTitleOffset(historg->GetXaxis()->GetTitleOffset());
    fHistogram->GetXaxis()->SetTitleFont(historg->GetXaxis()->GetTitleFont());

    fHistogram->GetYaxis()->SetTitle(historg->GetYaxis()->GetTitle());
    fHistogram->GetYaxis()->CenterTitle(historg->GetYaxis()->GetCenterTitle());
    fHistogram->GetYaxis()->RotateTitle(historg->GetYaxis()->GetRotateTitle());
    fHistogram->GetYaxis()->SetNoExponent(historg->GetYaxis()->GetNoExponent());
    fHistogram->GetYaxis()->SetNdivisions(historg->GetYaxis()->GetNdivisions());
    fHistogram->GetYaxis()->SetLabelFont(historg->GetYaxis()->GetLabelFont());
    fHistogram->GetYaxis()->SetLabelOffset(historg->GetYaxis()->GetLabelOffset());
    fHistogram->GetYaxis()->SetLabelSize(historg->GetYaxis()->GetLabelSize());
    fHistogram->GetYaxis()->SetTitleSize(historg->GetYaxis()->GetTitleSize());
    fHistogram->GetYaxis()->SetTitleOffset(historg->GetYaxis()->GetTitleOffset());
    fHistogram->GetYaxis()->SetTitleFont(historg->GetYaxis()->GetTitleFont());

    delete historg;
  }
  return fHistogram;

}
