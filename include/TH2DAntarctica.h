#ifndef TH2D_ANTARCTICA_H
#define TH2D_ANTARCTICA_H

#include "TProfile2D.h"
#include "TH2D.h"
#include "AntarcticaBackground.h"


/** 
 * Author: strutt@physics.ucla.edu
 * Date: September 2017
 * Description:
 * 
 * ANTARCTICA BACKGROUND ROOT-FRIENDLY HISTOGRAMS 
 * 
 * Easy to use! Just do something like:
 * 
 * TH2DAntarctica* h = new TH2DAntartica(name, title, nBinsX, nBinsY)
 * h->Fill(lon, lat, weight); // as you like, or, to test, h->FillRandomly(); 
 * h->Draw("colz");
 * 
 * The coordinate transformations are handled internally.
 * A pretty AntarcticaBackground is automatically drawn,
 * and via some laborious c-macro implemented magic,
 * can be maniplulated by the context menus of the histogram.
 * 
 * Sadly, actual implementation is quite ugly, 
 * and will be until ROOT does virtual inheritance.
 * I hope the simple usage will make up for that.
 * 
 * Happy histogramming.
 */




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
#define GF(var_type, internal_var_decl, FunctionName)                   \
  var_type Get##FunctionName() const                                    \
  {                                                                     \
    internal_var_decl f##FunctionName = getBackground()->Get##FunctionName(); \
    return f##FunctionName;                                             \
  }

#define SF(var_type, FunctionName)                         \
  void Set##FunctionName(var_type use##FunctionName)       \
  {                                                        \
    getBackground()->Set##FunctionName(use##FunctionName); \
  }










/** 
 * @class Class to draw pretty TProfile2Ds of Antarctica.
 * Handles coordinate transformations, so do Fill(longitude, latitude, value)
 * Puts a configurable AntarcticaBackground behind the image.
 */
class TProfile2DAntarctica : public TProfile2D {
 public:
  TProfile2DAntarctica(Int_t nx=-1, Int_t ny=-1);
  TProfile2DAntarctica(const char* name, const char* title, Int_t nx, Int_t ny);
  virtual ~TProfile2DAntarctica(){
    if(fAntarcticaBackground){
      delete fAntarcticaBackground;
      fAntarcticaBackground = NULL;
    }
  }

  virtual void Draw(Option_t* opt="");
  virtual Int_t Fill(Double_t lon, Double_t lat, Double_t val=1);

  void UnZoom(){getBackground()->UnZoom();} //*MENU

  GF(bool,bool,GrayScale)
  SF(bool,GrayScale)      //*TOGGLE *GETTER=GetGrayScale
  GF(bool,bool,ShowBases)
  SF(bool,ShowBases)      //*TOGGLE *GETTER=GetShowBases
  GF(bool,bool,ToolTip);
  SF(bool,ToolTip)        //*TOGGLE *GETTER=GetToolTip

  GF(int,,Coarseness)
  SF(int,Coarseness)      //*MENU *ARGS={useCoarseness=>fCoarseness}
  
  GF(bool,bool,Rampdem)
  SF(bool,Rampdem)        //*TOGGLE *GETTER=GetRampdem
  GF(bool,bool,Bed)
  SF(bool,Bed)            //*TOGGLE *GETTER=GetBed
  GF(bool,bool,Icemask)
  SF(bool,Icemask)        //*TOGGLE *GETTER=GetIcemask
  GF(bool,bool,Surface)
  SF(bool,Surface)        //*TOGGLE *GETTER=GetSurface
  GF(bool,bool,Thickness)
  SF(bool,Thickness)      //*TOGGLE *GETTER=GetThickness
  GF(bool,bool,Grid)
  SF(bool,Grid)           //*TOGGLE *GETTER=GetGrid
  
  void SetGridDivisions(Int_t deltaLon, Int_t deltaLat){
    getBackground()->SetGridDivisions(deltaLon, deltaLat);
  } // *MENU* *ARGS={deltaLat=>fDeltaLon, deltaLon=>fDeltaLat}
  
  virtual void ExecuteEvent(int event, int px, int py){
    getBackground()->ExecuteEvent(event, px, py);
    TProfile2D::ExecuteEvent(event, px, py);
  }
  void FillRandomly(Int_t nTimes = 5000);                               // *MENU
  void PrettifyPalette();                                               // *MENU
  void PrettifyBackgroundPalette(){getBackground()->PrettifyPalette();} // *MENU

 private:
  mutable AntarcticaBackground* fAntarcticaBackground; //! Don't persist
  mutable int                   fCoarseness;           //! Don't persist
  mutable int                   fDeltaLat;             //! Don't persist
  mutable int                   fDeltaLon;             //! Don't persist
  
  AntarcticaBackground* getBackground() const{
    if(!fAntarcticaBackground){
      fAntarcticaBackground = new AntarcticaBackground();
    }  
    return fAntarcticaBackground;
  }
  
  ClassDef(TProfile2DAntarctica, 1);

};









/** 
 * @class Class to draw pretty TH2Ds of Antarctica.
 * Handles coordinate transformations, so do Fill(longitude, latitude, weight)
 * Puts a configurable AntarcticaBackground behind the image.
 */

class TH2DAntarctica : public TH2D {
 public:
  TH2DAntarctica(Int_t nx=-1, Int_t ny=-1);
  TH2DAntarctica(const char* name, const char* title, Int_t nx, Int_t ny);
  virtual ~TH2DAntarctica(){
    if(fAntarcticaBackground){
      delete fAntarcticaBackground;
      fAntarcticaBackground = NULL;
    }
  }

  virtual void Draw(Option_t* opt="");
  virtual Int_t Fill(Double_t lon, Double_t lat, Double_t val=1);

  void UnZoom(){getBackground()->UnZoom();} //*MENU

  GF(bool,bool,GrayScale)
  SF(bool,GrayScale)      //*TOGGLE *GETTER=GetGrayScale
  GF(bool,bool,ShowBases)
  SF(bool,ShowBases)      //*TOGGLE *GETTER=GetShowBases
  GF(bool,bool,ToolTip);
  SF(bool,ToolTip)        //*TOGGLE *GETTER=GetToolTip

  GF(int,,Coarseness)
  SF(int,Coarseness)      //*MENU *ARGS={useCoarseness=>fCoarseness}
  
  GF(bool,bool,Rampdem)
  SF(bool,Rampdem)        //*TOGGLE *GETTER=GetRampdem
  GF(bool,bool,Bed)
  SF(bool,Bed)            //*TOGGLE *GETTER=GetBed
  GF(bool,bool,Icemask)
  SF(bool,Icemask)        //*TOGGLE *GETTER=GetIcemask
  GF(bool,bool,Surface)
  SF(bool,Surface)        //*TOGGLE *GETTER=GetSurface
  GF(bool,bool,Thickness)
  SF(bool,Thickness)      //*TOGGLE *GETTER=GetThickness
  GF(bool,bool,Grid)
  SF(bool,Grid)           //*TOGGLE *GETTER=GetGrid
  
  void SetGridDivisions(Int_t deltaLon, Int_t deltaLat){
    getBackground()->SetGridDivisions(deltaLon, deltaLat);
  } // *MENU* *ARGS={deltaLat=>fDeltaLon, deltaLon=>fDeltaLat}
  
  virtual void ExecuteEvent(int event, int px, int py){
    getBackground()->ExecuteEvent(event, px, py);
    TH2D::ExecuteEvent(event, px, py);
  }
  void FillRandomly(Int_t nTimes = 5000);                               // *MENU
  void PrettifyPalette();                                               // *MENU
  void PrettifyBackgroundPalette(){getBackground()->PrettifyPalette();} // *MENU

 private:
  mutable AntarcticaBackground* fAntarcticaBackground; //! Don't persist
  mutable int                   fCoarseness;           //! Don't persist
  mutable int                   fDeltaLat;             //! Don't persist
  mutable int                   fDeltaLon;             //! Don't persist
  
  AntarcticaBackground* getBackground() const{
    if(!fAntarcticaBackground){
      fAntarcticaBackground = new AntarcticaBackground();
    }  
    return fAntarcticaBackground;
  }
  
  ClassDef(TH2DAntarctica, 1);

};



#endif
