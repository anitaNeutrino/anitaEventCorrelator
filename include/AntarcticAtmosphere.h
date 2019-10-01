#ifndef _ANTARCTIC_ATMOSPHERE_H 
#define _ANTARCTIC_ATMOSPHERE_H 

/* Various routines for atmospheric modeling
 *
 * Some things depend on GeographicLib 
 *
 **/

#include "TGraph.h" 
#include "TF1.h" 
#include <set> 
#include <utility>


class Adu5Pat; 
namespace AntarcticAtmosphere
{

  /** Geoid model for converting from WGS84 to MSL and vice versa 
   * These require GeographicLib and the corresponding data files. See GeographicLib documentation for more details. 
   * */
  enum Geoid
  {
    EGM96_5, 
    EGM2008_1
  };
    
  /** Convert from mean sea level to WGS84 altitude. Requires GeographicLib */ 
  double MSLtoWGS84( double h, double lat, double lon, Geoid g = EGM96_5)  ; 

  /** Convert from WGS84 altitude to height above mean sea level. Requires GeographicLib*/ 
  double WGS84toMSL( const Adu5Pat * pat,  Geoid g = EGM96_5) ; 
  double WGS84toMSL( double lat, double lon, double alt,  Geoid g = EGM96_5) ; 


  /** Atmospheric Parameters. 
   * */ 
  struct Pars
  {
    double rho;  //kg / m^3 
    double P;    //mb 
    double T;    //Kelvin 
    double N;    //refractivity in units of (n-1) * 1e6 
  }; 

  enum Par
  {
    DENSITY,
    PRESSURE,
    TEMPERATURE,
    REFRACTIVITY
  };



  class AtmosphericModel 
  {
    public:   
      //h in m (msl), phi in radians (for position-dependent models. most models this has no effect) 
      virtual int computeAtmosphere(double h, Pars * p, double phi = 0) const= 0; 
      virtual double get(double h, Par p, double phi = 0) const; 
      TGraph * makeGraph(double hmin, double hmax, int nh, Par P, bool alt_on_x=true, double phi = 0) const; 
      virtual ~AtmosphericModel() { ; }
      virtual const char * name() const  { return "atmospheric model"; }
  }; 



  class InterpolatedAtmosphere : public AtmosphericModel
  {

    public: 
      InterpolatedAtmosphere() {;} 
      void addModel(const AtmosphericModel * m, double phi) { atms.insert(std::pair<double,const AtmosphericModel *>(phi,m));}
      virtual int computeAtmosphere(double h, Pars * p, double phi = 0) const; 
    private: 
      std::set<std::pair<double, const AtmosphericModel *> > atms; 
  }; 


  class SPRadiosonde : public AtmosphericModel
  {
    public: 

      /** This will download and use data from south pole radiosondes  */ 
      SPRadiosonde(int year, int month, int day, bool early = true); 

      virtual int computeAtmosphere(double h, Pars * p, double phi = 0) const; 
      
      double max_unextrapolated_height() const { return N.GetX()[N.GetN()-1]; }
      double min_height() const { return N.GetX()[0]; }
      bool ok() const { return loaded; } 
      virtual const char * name() const { return my_name.Data(); }

      const TGraph * raw_N() const { return &N ; }
 
    private: 
      bool loaded; 
      TGraph P; 
      TGraph T; 
      TGraph N; 
      TGraph rho; 
      TString my_name; 
      TF1 Nfit; 
      TF1 Pfit; 
      TF1 Tfit; 
      TF1 rhofit; 
  }; 






  class StandardUS : public AtmosphericModel
  {
    public: 
      StandardUS(double sea_level_T_kelvin = 265, double sea_level_P_mbar = 970)
      {
        sea_level_T = sea_level_T_kelvin; 
        sea_level_P = sea_level_P_mbar; 
      }
      virtual ~StandardUS() { ;} 

      virtual int computeAtmosphere(double h, Pars * p, double phi = 0) const; 
      virtual const char * name() const { return "Standard US"; }
    protected: 
      double sea_level_T; 
      double sea_level_P; 
  };

  /** Same as standard US atmosphere, but with refractivty = a exp(- b h) */ 
  class ExponentialRefractivity : public StandardUS
  {
    public: 
      ExponentialRefractivity(double a = 315, double b = 0.1361e-3) : StandardUS() {k_A = a ; k_B =b; my_name.Form("exponential refractivity, a = %g, b =%g", k_A, k_B); }

      virtual int computeAtmosphere(double h, Pars * p, double phi = 0) const; 
      virtual double get(double h, Par p, double phi = 0) const; 
      virtual ~ExponentialRefractivity() {; } 
      virtual const char * name() const { return my_name.Data(); }

    protected:
     double k_A, k_B; 
     TString my_name; 

  }; 

  const ExponentialRefractivity & ITURefractivity(); 




  class ArtificialInversion : public AtmosphericModel
  {
    public: 
      ArtificialInversion(const AtmosphericModel & base, double max_height, double max_ampl)  
        : m(base), hmax(max_height), Amax(max_ampl)  {; } 
      virtual int computeAtmosphere(double h, Pars * p, double phi = 0) const; 
      virtual double get(double h, Par p, double phi = 0) const; 
    private: 
      const AtmosphericModel & m;
      double hmax; 
      double Amax; 
      double correction(double h) const; 
  }; 
}



#endif
