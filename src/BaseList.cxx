#include "BaseList.h"
#include "AnitaVersion.h" 
#include "TFile.h" 
#include "TMath.h" 
#include <cmath>
#include "TTree.h" 
#include <unistd.h> 
#include "TROOT.h" 
#include "TKey.h" 


using namespace BaseList;


#ifndef DEG2RAD
#define DEG2RAD M_PI / 180
#endif


static void fillBases(std::vector<base> & baseList, int anita) 
{

  TString fname; 
  fname.Form("%s/share/anitaCalib/baseListA%d.root", getenv("ANITA_UTIL_INSTALL_DIR"), anita); 

  TFile fbase(fname.Data()); 

  if (!fbase.IsOpen())
  {
    fprintf(stderr,"Couldn't load base list for ANITA %d. Sorry :(\n", anita); 
    return;
  }

  //now load each tree 


  TIter iter(fbase.GetListOfKeys()); 
  TKey *k; 

  while ((k = (TKey *) iter()))
  {

    //only read in TTrees 
    TClass * cl = gROOT->GetClass(k->GetClassName()); 
    if (!cl->InheritsFrom("TTree")) continue; 

    TTree * t = (TTree*) k->ReadObj(); 
    TString source = t->GetName(); 

    std::string * str_name = 0; 
    double lon; 
    double lat; 
    double alt; 

    t->SetBranchAddress("name",&str_name); 
    t->SetBranchAddress("fullLat",&lat); 
    t->SetBranchAddress("fullLong",&lon); 
    t->SetBranchAddress("alt",&alt); 

    for (int i = 0; i < t->GetEntries(); i++) 
    {
      t->GetEntry(i); 
      baseList.push_back(base(TString(*str_name), source, lat,lon,alt)); 
    }
  }
}

static void fillPaths(std::vector<path> & pathList, int anita) 
{

  TString fname; 
  fname.Form("%s/share/anitaCalib/transientListRestrictedA%d.root", getenv("ANITA_UTIL_INSTALL_DIR"), anita); 

  //see if we have the restricted list

  if (access(fname.Data(),R_OK))
  {
    fprintf(stderr,"Couldn't find restricted list for ANITA %d (%s).  Will try to load unrestricted list. \n", anita, fname.Data()); 
    fname.Form("%s/share/anitaCalib/transientListUnrestrictedA%d.root", getenv("ANITA_UTIL_INSTALL_DIR"), anita); 
  }

  TFile fpath(fname.Data()); 

  if (!fpath.IsOpen())
  {
    fprintf(stderr,"Couldn't find unrestricted list for ANITA %d (%s).  Sorry :( \n", anita,  fname.Data()); 
    return; 
  }

  TIter iter(fpath.GetListOfKeys()); 
  TKey *k; 
  while ((k = (TKey *) iter()))
  {
  
    //only read in TTrees 
    TClass * cl = gROOT->GetClass(k->GetClassName()); 
    if (!cl->InheritsFrom("TTree")) continue; 

    TTree * t = (TTree*) k->ReadObj(); 

    TString source = t->GetName(); 
    std::string last_callsign = ""; 
    char callsign_buf[1024]; // 
    double lon; 
    double lat; 
    int alt = -1000; 
    int time; 

    // well,  I guess these trees are not as nicely normalized as the others. 


    t->SetBranchAddress("callSign",callsign_buf); 
    if (!t->GetBranch("fullLong")) //this tree has no position data. ignore it
    {
      continue ; 
    }


    t->SetBranchAddress("fullLong",&lon); 
    t->SetBranchAddress("fullLat",&lat); 
    t->SetBranchAddress("timeUTC",&time); 

    if (t->GetBranch("altitude")) // the traverse has no altitude data. Have no fear, we can fill it in ourselves. 
    {
      t->SetBranchAddress("altitude",&alt); 
    }

    std::vector<double> v_lat; 
    std::vector<double> v_lon; 
    std::vector<double> v_alt; 
    std::vector<unsigned> v_t; 

    for (int i = 0; i < t->GetEntries(); i++) 
    {
      t->GetEntry(i); 
      std::string callsign = callsign_buf; 

      if (last_callsign != callsign && v_lat.size()  )
      {
        pathList.push_back(path(TString(callsign.c_str()), source, (int) v_lat.size(), &v_lat[0], &v_lon[0], &v_alt[0], &v_t[0])); 
        v_lat.clear(); 
        v_lon.clear(); 
        v_alt.clear(); 
        v_t.clear(); 
      }

      last_callsign = callsign; 
      v_lat.push_back(lat); 
      v_lon.push_back(lon); 
      v_alt.push_back(alt); 
      v_t.push_back(time); 

    }

    //last iteration 
    if (v_lat.size()) 
    {
      pathList.push_back(path(TString(last_callsign.c_str()), source, v_lat.size(), &v_lat[0], &v_lon[0], &v_alt[0], &v_t[0])); 
    }
  }

}

// some annoying intermediate classes to be able to use magic statics 

static std::vector<base> no_bases; 
static std::vector<path> no_paths; 

struct baselist_impl 
{
  baselist_impl(int anita) 
  {
    fillBases(bases, anita); 
  }
  std::vector<base> bases; 

}; 

struct pathlist_impl 
{
  pathlist_impl(int anita) 
  {
    fillPaths(paths, anita); 
  }
  std::vector<path> paths; 

}; 



static std::vector<base> & bases()
{
  if (AnitaVersion::get() == 2) 
  {
    static baselist_impl bl(2); 
    return bl.bases; 
  }

  else if (AnitaVersion::get() == 3) 
  {
    static baselist_impl bl(3); 
    return bl.bases; 

  }
  else if (AnitaVersion::get() == 4) 
  {
    static baselist_impl bl(4); 
    return bl.bases; 
  }


  fprintf(stderr,"Don't have bases for %d\n", AnitaVersion::get()); 
  return no_bases; 
}

static std::vector<path> & paths()
{

  if (AnitaVersion::get() == 3) 
  {
    static pathlist_impl pl(3); 
    return pl.paths; 

  }
  else if (AnitaVersion::get() == 4) 
  {
    static pathlist_impl pl(4); 
    return pl.paths; 
  }

  fprintf(stderr,"Don't have paths for %d\n", AnitaVersion::get()); 
  return no_paths; 
}

const BaseList::base& BaseList::getBase(UInt_t index){

  index = index < bases().size() ? index : 0;
  return bases().at(index);
}

const BaseList::path& BaseList::getPath(UInt_t index){

  index = index < paths().size() ? index : 0;
  return paths().at(index);
}

const BaseList::abstract_base& BaseList::getAbstractBase(UInt_t index){

  if (index > bases().size() + paths().size()) index = 0; 
  return index < bases().size() ? (const BaseList::abstract_base &)  bases().at(index) : (const BaseList::abstract_base &) paths().at(index-bases().size()); 
}


size_t BaseList::getNumBases(){
  return bases().size();
}

size_t BaseList::getNumPaths() {
  return paths().size();
}

size_t BaseList::getNumAbstractBases(){
  return bases().size() + paths().size();
}



void BaseList::makeBaseList()
{
  makeEmptyBaseList(); 
  fillBases(bases(), AnitaVersion::get()); 
  fillPaths(paths(), AnitaVersion::get()); 
}


void BaseList::makeEmptyBaseList()
{
  bases().clear(); //DESTROY ALL THE BASES FOR SOME REASON  
}


BaseList::path::path(const TString & name, TString & source, 
                 int npoints, const double  * lat, const double * lon,
                 const double * alt,  const unsigned  * time)  


  : name(name) , dataSource(source), ts(time, time + npoints) 
{

  ps.reserve(npoints); 
  for (int i = 0; i < npoints; i++) 
  {
    ps.push_back(AntarcticCoord(AntarcticCoord::WGS84, lat[i], lon[i], alt[i]));
//    ps[i].to(AntarcticCoord::CARTESIAN);  //save as cartesian
  }


}

 
AntarcticCoord BaseList::path::getPosition(unsigned t) const {
 
  if (!isValid(t)) return AntarcticCoord(AntarcticCoord::WGS84, 90, 0, 0); // North pole is about as far as we can get! 

  //  Components to interpolate with.
  int l = TMath::BinarySearch(ts.size(), & ts [0], t); 
  int u = l + 1; 
  double low_frac = double(t - ts [l]) / double(ts [u] - ts [l]);  //  Lower fractional interpolative step.
  AntarcticCoord cl = ps [l];  //  Components in WGS84 coordinates.
  AntarcticCoord cu = ps [u];
  //  For later use; in case either altitude component isn't defined or negative. RampdemReader convention has longitude listed first, then latitude.
  double gndl = RampdemReader::SurfaceAboveGeoid(cl.y, cl.x, RampdemReader::surface);
  double gndu = RampdemReader::SurfaceAboveGeoid(cu.y, cu.x, RampdemReader::surface);
  double clz = (!cl.z || cl.z < 0) ? gndl : cl.z;
  double cuz = (!cu.z || cu.z < 0) ? gndu : cu.z;
  //  Great ellipse trajectories should correspond to the surface of the WGS84 geoid, so we should zero out out z-components.
  cl.z = 0, cu.z = 0;

  //  Cast vectors into Cartesian.
  cl.to(AntarcticCoord::CARTESIAN), cu.to(AntarcticCoord::CARTESIAN);

  //  Assuring a great ellipse trajectory between components, we will linearly interpolate polar gnomonic projection (X, Y, Z) between Cartesian points (x, y, z) ((X, Y, Z) = (1 / z) * (x, y, z)).
  //  See (https://www.uwgb.edu/dutchs/structge/sphproj.htm) and (http://mathworld.wolfram.com/StereographicProjection.html) for details.
  TVector3 g = low_frac / cu.z * cu.v() + (1 - low_frac) / cl.z * cl.v();
//  double rl = cl.v().Mag();
//  double ru = cu.v().Mag();
//  TVector3 g = low_frac * (ru / cu.z) * cu.v() + (1 - low_frac) * (rl / cl.z) * cl.v();

  //  Now to invert the transform, back to Cartesian ((x, y, z) = Z * (X, Y, Z), R = sqrt(X^2 + Y^2 + Z^2)).
  AntarcticCoord c = AntarcticCoord(g.z() * g);
//  double R = g.Mag();
//  AntarcticCoord c = AntarcticCoord((g.z() / R) * g);

  //  Return this Cartesian vector back in stereographic.
  if (clz == gndl && cuz == gndu) {

    c.to(AntarcticCoord::WGS84);
    double gnd = RampdemReader::SurfaceAboveGeoid(c.y, c.x, RampdemReader::surface);
    c.to(AntarcticCoord::STEREOGRAPHIC);
    c.z = gnd;
  } else {

    c.to(AntarcticCoord::STEREOGRAPHIC);
    c.z = low_frac * cuz + (1 - low_frac) * clz;  //  What we place as the stereographic z-component is actually the WGS84 component, altitude.
  }

  return c;

//  //  Interpolated components.
//  //  The following assumes spherical geometry, following a great circle trajectory. For our choice in Cartesian components, we must do things in terms of cotangents.
//  double cot_lat = low_frac / tan(DEG2RAD * cu.x) + (1 - low_frac) / tan(DEG2RAD * cl.x);
//  double lat = atan(1 / cot_lat) / DEG2RAD;
//  //  double lat = low_frac * cu.x + (1 - low_frac) * cl.x;
//  if (cu.y - cl.y < -180) cu.y += 360;  //  Accounting for longitude unwrapping, ensuring shorter longitude difference taken.
//  else if (cu.y - cl.y > 180) cu.y -= 360;
//  double lon = low_frac * cu.y + (1 - low_frac) * cl.y;
//  lon = fmod(lon + 180, 360) - 180;  //  Rewrapping longitude. Perhaps unneccessary if going to stereographically project anyway?
//  double alt = low_frac * cu.z + (1 - low_frac) * cl.z;
//  if (alt < 0) alt = RampdemReader::SurfaceAboveGeoid(lon, lat, RampdemReader::surface);  //  In case at least one of the input altitude components wasn't actually filled.
//
//  //  Construct the interpolated component vector, then return it stereographically projected.
//  AntarcticCoord c(AntarcticCoord::WGS84, lat, lon, alt);
//  c.to(AntarcticCoord::STEREOGRAPHIC);
//
//  return c;

//  AntarcticCoord cl = ps[l].as(AntarcticCoord::CARTESIAN); 
//  AntarcticCoord cu = ps[u].as(AntarcticCoord::CARTESIAN); 
//  double x =  low_frac * cl.x  + (1-low_frac) * cu.x; 
//  double y =  low_frac * cl.y  + (1-low_frac) * cu.y; 
//  double z =  low_frac * cl.z  + (1-low_frac) * cu.z; 
//
//  if (z < 0) {  //this means the altitude was not actually filled in (e.g for example a traverse), so we need to retrieve it ourselves... 
//
//    AntarcticCoord c(AntarcticCoord::CARTESIAN,x,y,0); 
//    c.to(AntarcticCoord::STEREOGRAPHIC);
//
//    c.z  = RampdemReader::SurfaceAboveGeoidEN(c.x,c.y, RampdemReader::surface); 
//    return c;
//  } else {
//    AntarcticCoord c(AntarcticCoord::WGS84, lat, lon, alt);
//    c.to(AntarcticCoord::STEREOGRAPHIC);
//    return c;
//  }
//  
//  //otherwise, we need to fix the altitude 
//
//  double alt = low_frac * ps[l].as(AntarcticCoord::WGS84).z + (1-low_frac)*ps[u].as(AntarcticCoord::WGS84).z; 
//
//  AntarcticCoord c(AntarcticCoord::WGS84, lat, lon, alt); 
//  AntarcticCoord c(AntarcticCoord::CARTESIAN,x,y,z); 
//  c.to(AntarcticCoord::STEREOGRAPHIC); //this is usually what we'll need
//  c.z = alt; //fix altitude 
//  return c; 
}



int BaseList::findBases(const char * query, std::vector<int> * matches, bool include_paths) 
{

  int first_found = -1; 
  for (unsigned i = 0; i < include_paths ? getNumAbstractBases() : getNumBases(); i++)
  {
    const abstract_base & a = getAbstractBase(i); 

    if (strcasestr(a.getName(), query))
    {
      if (first_found < 0) first_found = i; 
      if (matches)
      {
        matches->push_back(i); 
      }
      else break; 
    }
  }
  return first_found; 

} 
