#include "AntarcticAtmosphere.h" 
#include "TGraph.h" 
#include "TAxis.h" 
#include "Adu5Pat.h" 
#include "TObjArray.h" 
#include "TObjString.h" 
#include "TSystem.h" 
#include "math.h"
#include <algorithm> 

#ifdef USE_GEOGRAPHIC_LIB
#include "GeographicLib/Geoid.hpp" 
#endif 


const double REARTH = 6369; 
const double GMR = 34.163195; 


/**************************+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
****** US Atmosphere
 ***************************/


/*
 * It looks like US and INTL atmosphere are basically the same...
 *
 * much of this is from http://www.pdas.com/programs/atmos.f90 */ 
const int NTAB= 8; 
static double htab[] = {0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852  };
static double ttab[] = { 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946 };
static double ptab[] = { 1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6 }; 
static double gtab[] = {-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0 }; 

static int us_atmosphere(double alt_km, double * rel_rho,double * rel_P, double * rel_T)
{
  int i,j,k; 
  double h; 
  double tgrad, tbase; 
  double tlocal; 
  double deltah; 
  double sigma, delta, theta; 

  h = alt_km * REARTH / (alt_km + REARTH); 

  i = 0; 
  j = NTAB-1; 
 
  while (j > i+1) 
  {
    k = (i+j)/2; 

    if (h < htab[k])
    {
      j = k; 
    }
    else
    {
      i = k; 
    }
  }

  tgrad = gtab[i]; 
  tbase =ttab[i]; 
  deltah = h - htab[i]; 
  tlocal = tbase + tgrad * deltah; 

  theta = tlocal / ttab[0]; 

  if (tgrad == 0) 
  {
    delta = ptab[i] * exp(- GMR*deltah / tbase); 
  }
  else
  {
    delta = ptab[i] * pow(tbase / tlocal , GMR/tgrad); 
  }

  sigma = delta/theta; 

  *rel_rho = sigma; 
  *rel_P = delta; 
  *rel_T = theta; 

  return 0; 
}

int AntarcticAtmosphere::StandardUS::computeAtmosphere(double h, Pars *p, double phi) const 
{

  (void) phi; 

  int ret = us_atmosphere(h/1e3, &p->rho, &p->P, &p->T); 
  p->T *= sea_level_T; 
  p->P *= sea_level_P; 
  p->rho *= 1.225 * sea_level_T / sea_level_P; 
  p->N = 77.6 * p->P / p->T; 
  
  return ret; 
}

double AntarcticAtmosphere::ArtificialInversion::correction(double h) const
{
  if (h > hmax) return 1; 
  else return 1 - (hmax - h)/hmax * Amax; 

}


int AntarcticAtmosphere::ArtificialInversion::computeAtmosphere(double h, Pars *p, double phi) const 
{
  (void) phi; 
  int ret = m.computeAtmosphere(h,p); 
  p->T   *= correction(h);
  p->rho *= correction(h); 
  p->N   *= correction(h); 
  return ret; 
}

double AntarcticAtmosphere::ArtificialInversion::get(double h, Par p, double phi)  const
{
  (void) phi; 
  //short circuilt easy calculation 
  if (p == REFRACTIVITY) return m.get(h,p) * correction(h); 
  else if (p == TEMPERATURE) return m.get(h,p) * correction(h); 
  else if (p == DENSITY) return m.get(h,p) * correction(h); 
  else return m.get(h,p) ; 
}



int AntarcticAtmosphere::ExponentialRefractivity::computeAtmosphere(double h, Pars *p, double phi) const 
{
  (void) phi; 

  int ret = us_atmosphere(h/1e3, &p->rho, &p->P, &p->T); 
  p->T *= sea_level_T; 
  p->P *= sea_level_P; 
  p->rho *= 1.225 * sea_level_T / sea_level_P; 
  p->N = k_A * exp(-k_B * h); 
  return ret; 
}

double AntarcticAtmosphere::ExponentialRefractivity::get(double h, Par p, double phi)  const
{
  (void) phi; 
  //short circuilt easy calculation 
  if (p == REFRACTIVITY) return k_A * exp(-k_B* h) ; 
  Pars x; 
  computeAtmosphere(h,&x); 
  switch (p)
  {
    case DENSITY: 
      return x.rho;
    case PRESSURE: 
      return x.P;
    case TEMPERATURE: 
      return x.T;
    case REFRACTIVITY: 
      return x.N; 
    default: 
      return -9999; 
  }
}



#ifdef USE_GEOGRAPHIC_LIB

static const GeographicLib::Geoid & geoid2008_1() 
{
  static GeographicLib::Geoid  g("egm2008-1"); 
  return g; 
}

static GeographicLib::Geoid & geoid96_5() 
{
  static GeographicLib::Geoid  g("egm96-5"); 
  return g; 
}

static const GeographicLib::Geoid & getGeoid(AntarcticAtmosphere::Geoid g) 
{
  switch (g) 
  {
    case AntarcticAtmosphere::EGM2008_1: 
      return geoid2008_1(); 
    default:
      return geoid96_5(); 
  }
}

double AntarcticAtmosphere::MSLtoWGS84(double h, double lat, double lon, Geoid g)
{
  return  getGeoid(g).ConvertHeight ( lat, lon, h, GeographicLib::Geoid::GEOIDTOELLIPSOID); 
}

double AntarcticAtmosphere::WGS84toMSL(double lat, double lon, double alt, Geoid g)
{
  return getGeoid(g).ConvertHeight(lat,lon,alt, GeographicLib::Geoid::ELLIPSOIDTOGEOID); 
}

double AntarcticAtmosphere::WGS84toMSL(const Adu5Pat * pat, Geoid g)
{
  return  getGeoid(g).ConvertHeight ( pat->latitude, pat->longitude, pat->altitude, GeographicLib::Geoid::ELLIPSOIDTOGEOID); 
}

#else

double AntarcticAtmosphere::WGS84toMSL(const Adu5Pat * pat, Geoid g )
{
  static int nagged = 0; 

  if (nagged++ < 5) 
  {
    fprintf(stderr,"Need GeographicLib for MSL conversions\n"); 
  }
  return pat->altitude;

}


double AntarcticAtmosphere::WGS84toMSL(double lat, double lon, double alt, Geoid g) 
{

  static int nagged = 0; 

  if (nagged++ < 5) 
  {
    fprintf(stderr,"Need GeographicLib for MSL conversions\n"); 
  }
 
  return alt; 
}

double AntarcticAtmosphere::MSLtoWGS84(double h, double lat, double lon, Geoid g)
{
  static int nagged = 0; 

  if (nagged++ < 5) 
  {
    fprintf(stderr,"Need GeographicLib for MSL conversions\n"); 
  }
  return h;

}


#endif


TGraph * AntarcticAtmosphere::AtmosphericModel::makeGraph(double hmin, double hmax, int nh, Par p, bool alt_on_x, double phi) const
{ 

  TGraph * g = new TGraph(nh); 
  g->SetTitle(name()); 
  (alt_on_x ? g->GetYaxis() : g->GetXaxis())->SetTitle( p== DENSITY ? "Density (kg/m^3)" :
                           p== PRESSURE ? "Pressure (kPa) " : 
                           p== TEMPERATURE ? "Temperature (K)" : 
                           "(n  - 1) #times 10^{6}"); 

  (alt_on_x ? g->GetXaxis() : g->GetYaxis())->SetTitle("Height above MSL (m)"); 


  for (int i = 0; i < nh; i++)
  {
    double h = hmin + i*(hmax-hmin)/(nh-1); 

    double val = get(h,p,phi); 
    if (alt_on_x) 
      g->SetPoint(i, h,val); 
    else
      g->SetPoint(i, val,h); 
  }

  return g; 
}


double AntarcticAtmosphere::AtmosphericModel::get(double h, Par p, double phi)  const
{
  Pars x; 
  computeAtmosphere(h,&x,phi); 
  switch (p)
  {
    case DENSITY: 
      return x.rho;
    case PRESSURE: 
      return x.P;
    case TEMPERATURE: 
      return x.T;
    case REFRACTIVITY: 
      return x.N; 
    default: 
      return -9999; 
  }
}


const AntarcticAtmosphere::ExponentialRefractivity & AntarcticAtmosphere::ITURefractivity() 
{
 static ExponentialRefractivity er(315,0.1361e-3); 
 return er; 
}



//hopefully the format didn't change... 
static const char * sp_base = "ftp://amrc.ssec.wisc.edu/pub/southpole/radiosonde"; 


struct atm_meas
{
  double h; // m, msl
  double P; // hPa
  double T; // kelvin 
  double N; // refractivity (derived) 
  double rho; // density (derived) 

  bool operator< (const atm_meas & other) const
  {
    return h < other.h; 
  }
}; 



int AntarcticAtmosphere::SPRadiosonde::computeAtmosphere(double h, Pars *p, double phi) const
{
  (void) phi;
  p->rho = rho.Eval(h); 
  p->P = P.Eval(h); 
  p->T = T.Eval(h); 
  p->N = h > N.GetX()[N.GetN()-1] ? fit.Eval(h) : N.Eval(h); 
  return 0; 
}



AntarcticAtmosphere::SPRadiosonde::SPRadiosonde(int year, int mon, int day, bool early)
  : fit("Nfit","expo",10e3,50e3) 
{
  int late_time = 12; 
  TString cmd; 
  char * local_dir = getenv("RADIOSONDE_DIR"); 
  cmd.Form("%s %s/%d/%02d%02d%02d%02ddat.txt",local_dir ? "cat " : "curl", local_dir ? local_dir : sp_base,year,mon,day,year % 100, early? 0 :late_time); 
  TString data = gSystem->GetFromPipe(cmd.Data()); 
  loaded = false; 

  if (year < 2015) 
  {
    fprintf(stderr,"Haven't implemented parsing older radiosonde profiles yet...\n"); 
    return; 
  }

  my_name.Form("SP Radiosonde %d-%d-%d %s", year,mon,day, early ? " early " : "late"); 

  //split into lines 

  TObjArray  *lines = data.Tokenize("\n"); 
  lines->SetOwner(true); 

  // the south pole data format kindly changes from year to year... 


  std::vector<atm_meas> v; 


  //skip first 9 lines
  for (int i = 10; i < lines->GetEntries(); i++) 
  {
    const char * line = ((TObjString*)lines->At(i))->GetString().Data(); 

    int s,h,rh,dir; 
    float hPa,T,dwp,speed; 
    sscanf(line,"%d %d %g %g %g %d %g %d", &s,&h,&hPa,&T,&dwp,&rh,&speed,&dir);

    atm_meas m; 
    m.h = h; 
    m.P = hPa; 
    m.T= T + 273; 

    //I guess we assume it's ice? 

    double EF=1+1e-4 * (2.2 +m.P*(0.00382 + 6.4e-7 * (T*T))); 
    const double a = 6.115; 
    const double b = 23.036; 
    const double c = 279.82; 
    const double d = 333.7; 
    double e_s = EF * a * exp( (( b-T/d)*T)/(T+c)); 

    double e= e_s*rh/100.; 

    m.N = 77.6 * m.P / m.T - 5.6 * e / m.T + 3e5 * e / (m.T*m.T); 

    m.rho = (m.P) / (287.058 *m.T); //ignore wet for now since the air is plenty dry 
    v.push_back(m); 
     
  }

  std::sort(v.begin(), v.end()); 

  delete lines; 

  loaded = v.size(); 


  for (unsigned i = 0; i < v.size(); i++) 
  {
//    printf("%g %g %g %g, %g\n",v[i].h,v[i].P, v[i].T, v[i].N, v[i].rho); 
    P.SetPoint(i, v[i].h, v[i].P);
    T.SetPoint(i, v[i].h, v[i].T);
    N.SetPoint(i, v[i].h, v[i].N);
    rho.SetPoint(i, v[i].h, v[i].rho);
  }

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,7,0) 
  P.SetBit(TGraph::kIsSortedX); 
  T.SetBit(TGraph::kIsSortedX); 
  N.SetBit(TGraph::kIsSortedX); 
  rho.SetBit(TGraph::kIsSortedX); 
#endif

    N.Fit(&fit,"RQ"); 
}


int AntarcticAtmosphere::InterpolatedAtmosphere::computeAtmosphere(double h, Pars * p, double phi) const
{
  //find the two atmospheres to interpolate in between 

  if (atms.size()==0) return 1; 

  if (atms.size()==1) 
  {
    (*atms.begin()).second->computeAtmosphere(h,p,phi); 
    return 0; 
  }

  //otherwise let's find the bounds

  std::set<std::pair<double,const AtmosphericModel*> >::const_iterator it = atms.lower_bound(std::pair<double,const AtmosphericModel*>(phi,0)); 

  //this means our phi is smaller than the minimum
  if (it == atms.begin()) 
  {
    (*atms.begin()).second->computeAtmosphere(h,p,phi); 
    return 0; 
  }

  //this means our phi is bigger than maximum
  if (it == atms.end()) //our phi is greater than the last one 
  {

    (*atms.rbegin()).second->computeAtmosphere(h,p,phi); 
    return 0; 
  }


  //otherwise we want to interpolate between the two 
  double phi1 = (*it).first; 
  const AtmosphericModel * m1 =(*it).second;
  it--; 
  double phi0 = (*it).first; 
  const AtmosphericModel * m0 =(*it).second;

  double frac= phi1==phi0 ? 0.5: (phi-phi0)/(phi1-phi0); 

//  printf("%g | %g %g %g | %g\n", h, phi, phi0, phi1,frac); 
  Pars p0; 
  Pars p1; 

  m0->computeAtmosphere(h,&p0,phi);
  m1->computeAtmosphere(h,&p1,phi);

  p->rho = p0.rho * (1-frac) + frac*p1.rho; 
  p->P = p0.P * (1-frac) + frac*p1.P; 
  p->T = p0.T * (1-frac) + frac*p1.T; 
  p->N = p0.N * (1-frac) + frac*p1.N; 

  return 0; 

}


