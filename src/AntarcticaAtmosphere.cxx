#include "AntarcticAtmosphere.h" 
#include "TGraph.h" 
#include "TAxis.h" 
#include "Adu5Pat.h" 
#include "math.h"

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



int AntarcticAtmosphere::StandardUS::computeAtmosphere(double h, Pars *p) const 
{

  int ret = us_atmosphere(h/1e3, &p->rho, &p->P, &p->T); 
  p->T *= sea_level_T; 
  p->P *= sea_level_P; 
  p->rho *= 1.225 * sea_level_T / sea_level_P; 
  p->N = 77.6 * p->P / p->T; 
  
  return ret; 
}


#ifdef USE_GEOGRAPHIC_LIB

static const GeographicLib::Geoid & geoid2008_1() 
{
  static GeographicLib::Geoid  g("EGM2008-1"); 
  return g; 
}

static GeographicLib::Geoid & geoid96_5() 
{
  static GeographicLib::Geoid  g("EGM96-5"); 
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


TGraph * AntarcticAtmosphere::AtmosphericModel::makeGraph(double hmin, double hmax, int nh, Par p) const
{ 

  TGraph * g = new TGraph(nh); 
  g->GetXaxis()->SetTitle( p== DENSITY ? "Density (kg/m^3)" :
                           p== PRESSURE ? "Pressure (kPa) " : 
                           p== TEMPERATURE ? "Temperature (K)" : 
                           "(n  - 1) #times 10#{^6}"); 

  g->GetYaxis()->SetTitle("Height above MSL (m)"); 

  for (int i = 0; i < nh; i++)
  {
    double h = hmin + i*(hmax-hmin)/(nh-1); 

    g->SetPoint(i, get(h,p), h); 

  }

  return g; 
}


double AntarcticAtmosphere::AtmosphericModel::get(double h, Par p)  const
{
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



