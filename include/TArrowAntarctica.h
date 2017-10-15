#ifndef TARROW_ANTARCTICA
#define TARROW_ANTARCTICA

#include "TArrow.h"
#include "RampdemReader.h"


class TGraphAntarctica;

/**
 * @class TArrowAntarctica
 * @brief Draw an arrow between two lon/lat points
 * 
 * Handles conversion from lon/lat to Easting/Northing internally.
 * Currently doesn't have the fancy Draw behaviour that TGraphAntarctica does.
 * I suppose that could be added if required.
 */

class TArrowAntarctica : public TArrow {
public:

  /** 
   * Default constructor for ROOT, pretty useless
   */
  TArrowAntarctica () : TArrow () {setDefaultStyle();}

  /** 
   * Slightly useful constructor
   * 
   * @param lon1 is the longitude of the first point
   * @param lat1 is the latitude of the first point
   * @param lon2 is the longitude of the second point
   * @param lat2 is the latitude of the second point
   * @param arrowSize is the size of the arrow (default = 0.01)
   * @param option arrow draw option (default is "|>", which is a solid arrow head pointing towards the second point)
   */
  TArrowAntarctica(Double_t lon1, Double_t lat1, Double_t lon2, Double_t lat2, Float_t arrowSize=0.01, Option_t* option="|>")
    : TArrow (0, 0, 0, 0, arrowSize, option)
  {
    setDefaultStyle();
    SetPoint1(lon1, lat1);
    SetPoint2(lon2, lat2);
  }
  
  /** 
   * Perhaps an even more useful constructor.
   * Draws an arrow between two points of two TGraphAntarcticas (or the same one, if you pass it twice)
   * 
   * @param gr1 is the first TGraphAntarctica, defines the first point
   * @param gr2 is the second TGraphAntarctica, defines the second point
   * @param i1 is the point of the first TGraphAntarctica to use as the first point (default =0)
   * @param i2 is the point of the second TGraphAntarctica to use as the second point (default =0)
   * @param arrowSize is the size of the arrow (default = 0.01)
   * @param option arrow draw option (default is "|>", which is a solid arrow head pointing towards the second point)
   */
  TArrowAntarctica(TGraphAntarctica* gr1, TGraphAntarctica* gr2, Int_t i1 = 0, Int_t i2 = 0, Float_t arrowSize=0.01, Option_t* option="|>");
  

  /** 
   * Set the first arrow point
   * 
   * @param lon the longitude of the first point
   * @param lat the latitude of the first point
   */
  virtual void SetPoint1(Double_t lon, Double_t lat){
    RampdemReader::LonLatToEastingNorthing(lon, lat, fX1, fY1);
  }

  /** 
   * Set the second arrow point
   * 
   * @param lon the longitude of the second point
   * @param lat the latitude of the second point
   */  
  virtual void SetPoint2(Double_t lon, Double_t lat){
    RampdemReader::LonLatToEastingNorthing(lon, lat, fX2, fY2);
  }

private:
  /** 
   * Prettify by default
   */
  void setDefaultStyle(){
    SetLineWidth(2);
    EColor col = kRed;
    SetLineColor(col);
    SetFillColor(col);    
    SetOption("|>");
  }

  ClassDef(TArrowAntarctica, 1);
  
};


#endif //TARROW_ANTARCTICA
