#ifndef _ANITA_REFRACTION_MODEL_HH
#define _ANITA_REFRACTION_MODEL_HH

#include <vector> 
#include "AntarcticAtmosphere.h" 

class Adu5Pat; 
class AntarcticCoord; 
class TGraph; 
namespace Refraction
{

  /** Abstract interface for refraction models. 
   */ 
  class Model
  {

    public: 

     virtual ~Model() { ; } 
    /* Returns the elevation correction based on the source and payload position. This is the true elevation of the source
     * with respect to the apparent (true - apparent). 
     * Because I'm using the positive theta going down convention, this should be positive (the real theta is lower than the apparent). 
     * */ 
    virtual double getElevationCorrection(const Adu5Pat * pat, const AntarcticCoord * source, double * correction_at_source = 0) const = 0; 


  }; 



  /* A model that only depends on elevation, payload height and estimated source height */ 
  class PositionIndependentModel  : public Model
  {
    public: 
     virtual ~PositionIndependentModel() { ; } 
    /* Returns the elevation correction based on theta (positive down) , source height above MSL and payload height above MSL
     * 
     * */ 
    virtual double getElevationCorrection(double theta, double hSource, double hPayload, double * correction_at_source = 0) const = 0; 

    TGraph * makeGraph(double hSource, double hPayload, double minEl = 0); 

    /* Implemented in terms of the other prototpe */ 
    virtual double getElevationCorrection(const Adu5Pat * pat, const AntarcticCoord * source, double * correction_at_source = 0)const ; 
  }; 


  /** Peter's fits to his refraction code */ 
  class PGFit : public PositionIndependentModel
  {

    public: 
       virtual ~PGFit() { ; } 

      virtual double getElevationCorrection(double el, double hSource, double hPayload, double * correction_at_source = 0)const ; 
  };



  /** Raytracer  assuming Spherical Eart*/
  class RaytracerSpherical
  {
    public: 
      RaytracerSpherical(const AntarcticAtmosphere::AtmosphericModel *atm)
      {
        m = atm; 
        step_size = 100; 
      }

     struct Setup
     {
       double start_alt; 
       double end_alt; 
       double thrown_payload_angle; 
     };

     struct Result
     {
       double apparent_source_angle; 
       double actual_source_angle; 
       double actual_payload_angle; 
       double actual_distance; 
       double ray_distance; 
       double ray_time; 
       
     };

     int raytrace(const Setup * setup, Result * result); 

     double step_size; 

     int nSteps() const { return (int) last_path_x.size(); }
     const double * lastY() { return &last_path_y[0]; } 
     const double * lastX() { return &last_path_x[0]; } 

     void draw(); 

    private: 
     const AntarcticAtmosphere::AtmosphericModel *m; 
     std::vector<double> last_path_x; 
     std::vector<double> last_path_y; 
  }; 


}

#endif
