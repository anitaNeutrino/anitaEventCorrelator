#ifndef TPROFILE2D_ANTARCTICA_H
#define TPROFILE2D_ANTARCTICA_H

#include "TProfile2D.h"
#include "AntarcticaBackground.h"



/** 
 * Silly pre-processor macros to wrap calls to the background.
 * What's the point?
 * To give you the lovely interactive right-click context,
 * menus in ROOT without doing too much typing.
 * 
 * To use the interative features of the background, do
 * 
 * GF(bool,GrayScale);
 * SF(bool,GrayScale) // *TOGGLE *GETTER=GetGrayScale
 * 
 * It would be better to do this in a single macro, 
 * but the preprocessor can't do newlines and rootcling
 * attaches the TOGGLE/MENU nonsense to each function defined
 *  on the same line, hence the need for two macros taking 
 * identical arguments, and the comment string after SF.
 * @param var_type is the type of the getter/setter function
 * @param FunctionName is the name of the function in AntarcticaBackground
 */
#define GF(var_type, FunctionName)               \
  var_type Get##FunctionName() const             \
  {                                              \
    return getBackground()->Get##FunctionName(); \
  }

#define SF(var_type, FunctionName)                         \
  void Set##FunctionName(var_type use##FunctionName)       \
  {                                                        \
    getBackground()->Set##FunctionName(use##FunctionName); \
  }






class TProfile2DAntarctica : public TProfile2D {
 public:
  TProfile2DAntarctica(Int_t nx=-1, Int_t ny=-1);
  virtual ~TProfile2DAntarctica();

  virtual void Draw(Option_t* opt="");

  virtual Int_t Fill(Double_t lon, Double_t lat, Double_t weight=1);

  GF(bool,GrayScale)
  SF(bool,GrayScale)      //*TOGGLE *GETTER=GetGrayScale
  GF(bool,ShowBases)
  SF(bool,ShowBases)      //*TOGGLE *GETTER=GetShowBases
  GF(bool,ToolTip);
  SF(bool,ToolTip)        //*TOGGLE *GETTER=GetToolTip
  
  GF(bool,Rampdem)
  SF(bool,Rampdem)        //*TOGGLE *GETTER=GetRampdem
  // GF(bool,Bed)
  // SF(bool,Bed)            //*TOGGLE *GETTER=SetBed
  // GF(bool,Icemask)
  // SF(bool,Icemask)        //*TOGGLE *GETTER=SetIcemask
  // GF(bool,Surface)
  // SF(bool,Surface) //*TOGGLE *GETTER=SetSurface
  // GF(bool,Thickness)
  // SF(bool,Thickness) //*TOGGLE *GETTER=SetThickness
  // GF(bool,Grid)
  // SF(bool,Grid) //*TOGGLE *GETTER=SetGrid
  
  // void SetGridDivisions(Int_t deltaLon, Int_t deltaLat){    
  //   getBackground()->SetGridDivisions(deltaLon, deltaLat);
  // } // *MENU* *ARGS={deltaLat=>fDeltaLon, deltaLon=>fDeltaLat}  
  
  virtual void ExecuteEvent(int event, int px, int py);  
  void FillRandom(Int_t nTimes = 5000);

 private:
  mutable AntarcticaBackground* fAntarcticaBackground; //! Don't persist
  AntarcticaBackground* getBackground() const;
  
  ClassDef(TProfile2DAntarctica, 1);

};


#endif
