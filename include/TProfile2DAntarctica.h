#ifndef TPROFILE2D_ANTARCTICA_H
#define TPROFILE2D_ANTARCTICA_H

#include "TProfile2D.h"

class AntarcticaBackground;

class TProfile2DAntarctica : public TProfile2D {
 public:
  TProfile2DAntarctica(Int_t nx=-1, Int_t ny=-1);
  virtual ~TProfile2DAntarctica();

  virtual void Draw(Option_t* opt="");

  virtual Int_t Fill(Double_t lon, Double_t lat, Double_t weight=1);

  void FillRandom(Int_t nTimes = 5000);

  AntarcticaBackground* getBackground(); // generate the background  
 private:
  AntarcticaBackground* fAntarcticaBackground; //! Does not persist in ROOT!
  
  ClassDef(TProfile2DAntarctica, 1);


};


#endif
