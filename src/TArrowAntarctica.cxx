#include "TArrowAntarctica.h"
#include "TGraphAntarctica.h"

ClassImp(TArrowAntarctica);

TArrowAntarctica::TArrowAntarctica(TGraphAntarctica* gr1, TGraphAntarctica* gr2, Int_t i1, Int_t i2, Float_t arrowSize, Option_t* option)
  : TArrow (0, 0, 0, 0, arrowSize, option)
{
  setDefaultStyle();
  if(gr1 && i1 >= 0 && i1 < gr1->GetN() ){
    fX1 = gr1->GetX()[i1];
    fY1 = gr1->GetY()[i1];
  }
  if(gr2 && i2 >= 0 && i2 < gr2->GetN() ){
    fX2 = gr2->GetX()[i2];
    fY2 = gr2->GetY()[i2];
  }
}
