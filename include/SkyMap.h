#ifndef SKYMAP_HH
#define SKYMAP_HH

#include <vector> 
#include "TH2.h" 
#include "TPaletteAxis.h" 
#include "TMarker.h" 
#include "TGraph.h" 
#include "TText.h" 
#include "TLatex.h" 

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

    SkyMap(const SkyMap & other); 

    virtual void Paint(Option_t * option = ""); 
    void setTitle(const char * title) { sky_background.SetTitle(title) ; } 
    void setReverseX(bool r) { reverseX = r; } 
    void setBackground(const TH2 * background, bool ra_in_hours=false); 
    void addMarker(const TMarker * marker); 
    void addGraph(const TGraph * graph, bool ra_in_hours=false); 
    void addText(const TLatex * txt); 
    void clearMarkers() { sky_markers.clear() ; }
    void setMinimum(double min) { sky_background.SetMinimum(min); }
    void setMaximum(double max) { sky_background.SetMaximum(max); }
    TH2 & getBackground() { return sky_background ; } 
    void clearGraphs() { sky_graphs.clear() ; } 
    void ExecuteEvent(Int_t event, Int_t px, Int_t py) { sky_background.ExecuteEvent(event,px,py); } 

    void setRange(double min_ra, double max_ra, double min_dec, double max_dec); 
    void unsetRange() { range_set = false; }
    virtual ~SkyMap(); 

    static void toMollweide(double lon, double lat, double &x, double &y,double lon_0 = 0,  bool reverse_x = true); 
    static void fromMollweide(double x, double y, double & lon, double &lat, double lon_0 = 0, bool reverse_x = true); 

  private: 
    
    bool reverseX; 
    double lon_0; 
    bool range_set; 
    double range[4]; 
    TH2D sky_background; 
    std::vector<TMarker> sky_markers; 
    std::vector<TGraph> sky_graphs; 
    std::vector<TLatex> sky_text; 
    TText left; 
    TText right; 
    TPaletteAxis * zaxis; 
    ClassDef(SkyMap,1); 

};  



#endif

