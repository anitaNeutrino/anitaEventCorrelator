#ifndef ANTARCTICA_GEOM_H
#define ANTARCTICA_GEOM_H



/** Some things to help divide up Antarctica 
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 * */ 

#include "RampdemReader.h" 
#include <vector> 
#include "TVector3.h" 

class Adu5Pat; 



/** Slightly smarter Antarctic Coordinates.  */ 
class AntarcticCoord
{
  public: 

    enum CoordType
    {
      WGS84,  // (lat lon alt)  [deg, deg, m] 
      STEREOGRAPHIC,  // (Easting Northing alt) [m,m,m]  
      CARTESIAN  // (x y z)  [m,m,m] 
    }; 


    AntarcticCoord(CoordType coord_type = CARTESIAN, double xval = 0, double yval = 0, double zval = 0)
    {
      set (coord_type,xval,yval,zval) ; 
    }

    AntarcticCoord(const TVector3 & v) 
    {
      set(CARTESIAN, v.x(),v.y(),v.z()); 
    }

    TVector3 v() const
    {
      AntarcticCoord c = as(CARTESIAN); 
      return TVector3(c.x,c.y,c.z); 
    }

    virtual ~AntarcticCoord() { ; } 

    void set(CoordType coord_type , double xval , double yval , double zval )
    {
      type = coord_type; 
      x= xval; 
      y= yval; 
      z= zval; 
    }



    /* Get a copy of this AntarcticCoord with the chosen coord system */  
    const AntarcticCoord as(CoordType new_type) const { AntarcticCoord c(*this); c.to(new_type); return c; }

    /* Convert this coordinate to the type you want */ 
    void to(CoordType new_type)  { if (new_type != type) {convert(new_type);} }

    CoordType currentType() const { return type; } 
    double x; 
    double y; 
    double z; 

    void asString(TString * str) const; 
  private:
    CoordType type; 
    void convert(CoordType new_type); 

    ClassDef(AntarcticCoord,1); 

}; 


/* A segmentation scheme divides Antarctica into chunks (segments) 
 *
 * Each piece has a unique index, which can be accessed from its coordinates.
 * The center of each segment and samples from it may be queried from the
 * index. This will need altitude information, the source of which can be
 * changed with setRampdemDataset(). 
 *
 */ 

class AntarcticSegmentationScheme {

  public: 

  /** Factory for AntarcticSegmentationScheme . Will create one based on a
   * string key. 
   *
   * Currently implemented:
   * stereographic_grid_nx[I]_ny[I,=nx]_maxE[I,=3330000]_maxN[I,=3330000]  nx
   * by ny bin segmentation of antarctica 
   *
     
   *
   * Planned: 
   *  healpix
   *
   */ 
    static AntarcticSegmentationScheme * factory(const char * key); 

    virtual ~AntarcticSegmentationScheme() { ; } 

    /** Gets the index of the segment at this position. -1 if not in segmentation scheme.  */ 
    virtual int getSegmentIndex(const AntarcticCoord & coord) const= 0; 

    virtual int NSegments() const = 0;  

    /* Gets the center of the segment at this position.
     * 
     * If fillalt, will also fill altitude using RampdemReader (which will be somewhat slower) , otherwise sets it to equivalent of 0 
     *
     * */
    virtual void getSegmentCenter(int idx, AntarcticCoord * fillme, bool fillalt = true) const = 0; 

    /** convenience method */ 
    virtual AntarcticCoord getSegmentCenter(int idx, bool fillalt = true) const { AntarcticCoord c; getSegmentCenter(idx,&c,fillalt); return c; } 

    /** Return the number of neighbhors. 
     * if neighbhors not null, appends to the vector the  neighbors. 
     *
     */ 
    virtual int getNeighbors(int segment, std::vector<int> * neighbors = NULL) const = 0; 

    /** Sample N coordinates within this segment. Will allocate N coords fillus is null. If random is true, random sampling, otherwise some uniformish type 
     *
     * If fillalt, will also fill altitude using RampdemReader (which will be somewhat slower) , otherwise sets it to equivalent of 0 
     * */ 
    virtual AntarcticCoord * sampleSegment(int idx, int N, AntarcticCoord * fillus = 0, bool random = true, bool fillalt = true ) const = 0; 

    /** Draw the segmentation scheme in Stereographic Coords. Default implementation samples each segment 64 times  and draws as a TGraph2D with z = segment numbe
     *
     *
     * if data is 0, draws the index of each segment, otherwise draws the value in data[i] for each segment.
     * */
    virtual void Draw(const char * opt = "colz", const double * data = 0 ) const; 
    virtual void DrawI(const char * opt = "colz", const int * data = 0) const; 



    virtual void asString(TString * str) const = 0; 
    virtual void setRampdemDataset( RampdemReader::dataSet d)  { dataset = d; } 

  protected: 
    AntarcticSegmentationScheme() :  dataset (RampdemReader::rampdem) { ;} 
    RampdemReader::dataSet dataset; 

    ClassDef(AntarcticSegmentationScheme,1); 

}; 




class StereographicGrid : public AntarcticSegmentationScheme
{

  public:
    StereographicGrid(int nx = 100, int ny = 100,  double max_E = 3330000, double max_N = 3330000); 
    virtual ~StereographicGrid()  {; } 
    virtual int getSegmentIndex(const AntarcticCoord & coord) const; 
    virtual int getNeighbors(int segment, std::vector<int> * neighbors = NULL) const; 
    virtual int NSegments() const { return nx * ny; } 
    virtual void getSegmentCenter(int idx, AntarcticCoord * fill, bool fillalt=true) const; 
    virtual AntarcticCoord * sampleSegment(int idx, int N, AntarcticCoord * fillus = 0, bool random = true, bool fillalt = true) const; 
    virtual void Draw(const char * opt = "colz", const double * data = 0) const; 

    virtual void asString(TString * str) const; 
  private:

    int nx; 
    int ny; 
    double max_E; 
    double max_N; 
    double dx; 
    double dy; 

    ClassDef(StereographicGrid,1); 


}; 


/** 
 * Class to hold angular positions of source / payload in each other's frame
 *
 */ 
class PayloadParameters
{
  public: 
    PayloadParameters(const Adu5Pat * gps, const AntarcticCoord & source_pos); 
    double source_phi;  //the phi of the source, in payload coordinates (degrees). 
    double source_theta; //the theta of the source, in payload coordinates (degrees) such that theta > 0 is coming from below
    double payload_el;  // the elevation of the payload from the source (deg) . positive is UP
    double payload_az; // the azimuth of the payload from the source (deg). 
    double distance; //distance between source and payload
    AntarcticCoord payload; 
    AntarcticCoord source; 
  private: 

    PayloadParameters();  
    ClassDefNV (PayloadParameters,1); 
};


/*
class HealPixSegmentation : public AntarcticSegmentationScheme
{

  public: 
    HealPixSegmentation(int npix, double max_lat = 60); 
    virtual ~HealPixSegmentation(); 
    virtual int getSegmentIndex(const AntarcticCoord & coord) const; 
    virtual int NSegments() const { return nx * ny; } 
    virtual void getSegmentCenter(int idx, AntarcticCoord * fill, bool fillalt=true) const; 
    virtual AntarcticCoord * sampleSegment(int idx, int N, AntarcticCoord * fillus = 0, bool random = true, bool fillalt = true) const; 
    virtual void asString(TString * str) const; 
  private:

}; 

*/

#endif
