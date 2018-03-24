#ifndef SKYMAP_HH
#define SKYMAP_HH

#include <vector> 
#include "TH2.h" 
#include "TPaletteAxis.h" 
#include "TMarker.h" 
#include "TGraph.h" 
#include "TText.h" 

/* My attempt at making a sky map. It supports: 
 *   
 *   A background TH2 
 *   Any number of TGraphs ("contours") and TMarkers
 *
 *   All input units are in degrees. 
 *
 */ 


class SkyMap : public TObject
{

  public: 
    SkyMap(double lon_0 = 180, int nbinsx = 1000, int nbinsy = 1000, const TH2 * background = 0, const std::vector<const TMarker *> *  markers = 0, 
                                   const std::vector<const TGraph*>  * graphs = 0); 

    virtual void Paint(Option_t * option = ""); 
    void setTitle(const char * title) { sky_background.SetTitle(title) ; } 
    void setReverseX(bool r) { reverseX = r; } 
    void setBackground(const TH2 * background); 
    void addMarker(const TMarker * marker); 
    void addGraph(const TGraph * graph); 
    void clearMarkers() { sky_markers.clear() ; }
    void clearGraphs() { sky_graphs.clear() ; } 

    virtual ~SkyMap(); 

    static void toMollweide(double lon, double lat, double &x, double &y,double lon_0 = 0,  bool reverse_x = true); 
    static void fromMollweide(double x, double y, double & lon, double &lat, double lon_0 = 0, bool reverse_x = true); 

  private: 
    
    bool reverseX; 
    double lon_0; 
    TH2D sky_background; 
    std::vector<TMarker> sky_markers; 
    std::vector<TGraph> sky_graphs; 
    TText left; 
    TText right; 
    TPaletteAxis * zaxis; 
    ClassDef(SkyMap,1); 

};  



#endif

