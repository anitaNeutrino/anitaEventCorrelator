#include "AntarcticaGeometry.h"
#include "RefractionModel.h"

class Adu5Pat;

/** 
 * Class to hold angular positions of source / payload in each other's frame
 *
 */ 
class PayloadParameters
{
  public: 
    PayloadParameters(const Adu5Pat * pat,  const AntarcticCoord & source_pos, const Refraction::Model * refraction =0); 
    PayloadParameters () { ; }
    PayloadParameters(const PayloadParameters & other); 

    double source_phi;  //the phi of the source, in payload coordinates (degrees). 
    double source_theta; //the theta of the source, in payload coordinates (degrees) such that theta > 0 is coming from below
    double payload_el;  // the elevation of the payload from the source (deg) . positive is UP
    double payload_az; // the azimuth of the payload from the source (deg). 
    double distance; //distance between source and payload

    /* Checks for collision with the ground */
    bool checkForCollision(double dx = 100, AntarcticCoord * where = 0, AntarcticCoord * where_exit = 0, RampdemReader::dataSet d = RampdemReader::rampdem, double grace = 20, bool reverse = false) const;

    AntarcticCoord payload; 
    AntarcticCoord source; 

    /** This is my version of the trace back to continent functions. 
     * I'm not sure it does exactly the right thing because of the altitude change, but it
     * follows the geodesic in the phi direction until the payload coordinates match
     * what they're supposed to. 
     *
     * Returns 1 on success, 0 if over horizon (but fills in payload position with the point at the horizon), -1 if doesn't converge to tolerance (but fills in closest) 
     *
     */

    static int findSourceOnContinent(double theta, double phi,
                          const Adu5Pat * gps, PayloadParameters * fillme, 
                          const Refraction::Model * m = 0, 
                          double collision_check_dx = 0, // 0 to not cehck 
                          double min_dx = 5,
                          double tol = 5e-6, 
                          double min_el  = 0,
                          RampdemReader::dataSet d= RampdemReader::rampdem); 

  private: 

    ClassDefNV (PayloadParameters,1); 
};



