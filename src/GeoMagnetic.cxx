#include "GeoMagnetic.h"

#include <iostream>
#include <fstream>
#include "TString.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TObjString.h"
#include <map>
#include "RampdemReader.h"
#include "AntarcticaBackground.h"
#include "TF1.h"

#include "AnitaGeomTool.h"
#include "TDatime.h"

// --------------------------------------------------------------------------------------------------------------------------------------
// Silly globals, but I want to hide the internal complexities inside the namespace
// --------------------------------------------------------------------------------------------------------------------------------------

bool doneInit = false; // Tells you whether we've read in the data and precalculated the factorials
const int numPoly = 14; // actually there are only 12 polynomial coeffients, but I'm going to start counting from one for simplicity
std::vector<double> factorials(2*numPoly, 0);
std::map<int, std::vector<double> > g_vs_time; // Gauss coefficients (needed to calc potential)
std::map<int, std::vector<double> > h_vs_time; // Gauss coefficients (needed to calc potential)
const double earth_radius = 6371.2e3; // earth radius in meters for the magnetic model
TF1* fAssocLegendre[numPoly][numPoly] = {{NULL}}; // Associated Legendre polynomials



// for differentiating the potential
const double dr = 1;
const double dTheta = 0.01*TMath::DegToRad();  
const double dPhi = 0.01*TMath::DegToRad();

typedef enum {
  r_hat,
  theta_hat,
  phi_hat
} unitVec;



// --------------------------------------------------------------------------------------------------------------------------------------
// Utility functions, for initialisation and internal calculation
// --------------------------------------------------------------------------------------------------------------------------------------

TVector3 makeUnitVectorSpherical(double theta, double phi, unitVec dir){
  TVector3 unit;
  switch(dir){
    case r_hat:
      unit.SetXYZ(TMath::Sin(theta)*TMath::Cos(phi),
                  TMath::Sin(theta)*TMath::Sin(phi),
                  TMath::Cos(theta));
    break;    
    
    case theta_hat:
      unit.SetXYZ(TMath::Cos(theta)*TMath::Cos(phi),
                  TMath::Cos(theta)*TMath::Sin(phi),
                  -TMath::Sin(theta));
    break;
    case phi_hat:
      unit.SetXYZ(-TMath::Sin(phi),
                  TMath::Cos(phi),
                  0);
      break;
    default:
      std::cerr << "Error in " << __FILE__ << ", unknown unit vector requested..." << std::endl;
      break;
  }
  return unit;
}

/** 
 * Maps the polynomials degree(n) and order (n) to an index for the vectors inside the g_vs_time, g_vs_time maps
 * 
 * @param n is the degree
 * @param m is the order
 * 
 * @return the vector index
 */
inline double getIndex(int n, int m){
  return n*numPoly + m;
}

/** 
 * Reads in the Gauss coefficients for the associated Legendre polynomials
 * 
 */
void getGaussCoefficients(){  

  const char* anitaUtilInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  if(!anitaUtilInstallDir){
    std::cerr << "Warning in " << __FILE__ << ", ANITA_UTIL_INSTALL_DIR not set" << std::endl;
  }
  TString fileName = TString::Format("%s/share/anitaCalib/igrf12coeffs.txt", anitaUtilInstallDir);
  std::ifstream coeffs(fileName);

  if(!coeffs.good()){
    std::cerr << "Error in " << __FILE__ << " unable to find " << fileName.Data() << std::endl;
  }

  const int expectedTokens = 28;  
  std::vector<std::vector<TString> > igrfDataTableStrings;
  while(!coeffs.eof()){
    std::string line;
    std::getline(coeffs, line);
    
    if(line[0] != '#'){ // ignore comment lines at start

      TString s =  line.c_str();

      // remove windows style newline prefix (for \r\n vs. \n)
      s.ReplaceAll("\r", "");

      TObjArray* tokens = s.Tokenize(' ');
      int numTokens = tokens->GetEntries();

      if(numTokens > 0){ // last line has nothing in
        igrfDataTableStrings.push_back(std::vector<TString>());
        
        if(numTokens!=expectedTokens){
          std::cerr << "Warning parsing " << fileName << ", found " << numTokens << " tokens in line..." << std::endl;
          std::cerr << line << std::endl;
        }

        for(int i=0; i < numTokens; i++){
          TString thisS = ((TObjString*)tokens->At(i))->GetString();          
          igrfDataTableStrings.back().push_back(thisS);
        }
      }
    }
  }

  
  // now extract data from strings...
  const int numRows = igrfDataTableStrings.size();

  // there are two rows of column titles so start at third row
  // std::vector<TString>& headerRow1 = igrfDataTableStrings.at(0);
  std::vector<TString>& headerRow2 = igrfDataTableStrings.at(1);

  // extract years from the header
  std::map<unsigned, int> colToYear;
  for(unsigned col=3; col < headerRow2.size() - 1; col++){
    int year = atoi(headerRow2[col].Data());
    colToYear[col] = year;
    // std::cout << col << "\t" << year << std::endl;
  }

  //  loop over rows and extract Gauss coefficents
  for(int row=2; row < numRows; row++){
    std::vector<TString>& thisRow = igrfDataTableStrings.at(row);

    // first three rows tell you whether it's g or h, and the m and n Legendre degree and order
    const TString& g_or_h = thisRow[0];
    const int n = atoi(thisRow[1].Data());
    const int m = atoi(thisRow[2].Data());
    
    int index = getIndex(n, m);
    
    // std::cout << g_or_h << "\t" << m << "\t" << n << "\t" << index << std::endl;
    
    for(unsigned col=3; col < thisRow.size() - 1; col++){
      int year = colToYear[col];

      std::map<int, std::vector<double> >::iterator it;      

      // put data into Gauss coeffient maps
      if(g_or_h == "g"){

        it = g_vs_time.find(year);
        if(it!=g_vs_time.end()){
          it->second.at(index) = atof(thisRow.at(col).Data());
        }
        else{
          g_vs_time[year] = std::vector<double>(numPoly*numPoly, 0);
          g_vs_time[year].at(index) = atof(thisRow.at(col).Data());
        }
      }
      
      else if (g_or_h == "h"){
        it = h_vs_time.find(year);
        if(it!=h_vs_time.end()){
          it->second.at(index) = atof(thisRow.at(col).Data());
        }
        else{
          h_vs_time[year] = std::vector<double>(numPoly*numPoly, 0);
          h_vs_time[year].at(index) = atof(thisRow.at(col).Data());
        }
      }
      else{
        std::cerr << "Warning! Unknown coefficient " << g_or_h << "[" << n << "][" << m << "]" << std::endl;
      }
    }
  }
}







/** 
 * Do all the precalculation and initialisation needed to make the namespace useable
 *
 * This involves reading in the coefficients, precalculating the factorials
 * and making the TF1 associated legendre polynomials
 */
void init(){

  if(!doneInit){
    getGaussCoefficients();
    for(int i=0 ; i< 2*numPoly; i++){
      factorials.at(i) = TMath::Factorial(i);
      // std::cout << i << "! = " << factorials.at(i) << std::endl;
    }
    for(int n=1; n < numPoly; n++){
      for(int m=0; m <=n; m++){
        if(!fAssocLegendre[n][m]){
          TString name = TString::Format("fAssocLegendre_%d_%d", n,  m);    
          TString formula = TString::Format("ROOT::Math::assoc_legendre(%d, %d, x)", n, m);
          fAssocLegendre[n][m] = new TF1(name, formula, -1, 1);
          // std::cout << m << "\t" << n << std::endl;
        }
      }
    }  
    doneInit = true;
  }
}


/** 
 * Look up precalculated factorial, prints warning if outside precalculated range
 *
 * ROOT may do this internally, who knows? Perhaps this is premature optimisation...
 * 0! factorial is stored in factorials[0], and so on up to 27! in factorials[27].
 *
 * @param i is the number you wish you find the factorial for, must be less than 28
 * 
 * @return i factorial (or -1 if i < 0, i > 27)
 */
double getFactorial(int i){
  if(i >=  2*numPoly){
    std::cerr << "Too high factorial requested!" << std::endl;
    return -1;
  }    
  return factorials[i];
}




/** 
 * Convert unixTime to fractional year with waaayyy too much precision
 * I think this should correctly handle leap years and other anomolies
 * 
 * @param unixTime is the seconds since 1970
 * 
 * @return year as a decimal quantity
 */
double unixTimeToFractionalYear(UInt_t unixTime){
  TDatime t2(unixTime);
  int thisYear = t2.GetYear();
  TDatime t1(thisYear, 0, 0, 0, 0, 0);
  UInt_t unixTimeYearStart = t1.Convert();
  TDatime t3(thisYear+1, 0, 0, 0, 0, 0);
  UInt_t unixTimeNextYear = t3.Convert();
  double year = thisYear;
  year += double(unixTime - unixTimeYearStart)/double(unixTimeNextYear - unixTimeYearStart);
  // std::cout << unixTime << "\t" << thisYear << "\t" << unixTimeYearStart << "\t" << unixTimeNextYear << std::endl;
  return year;
}




/** 
 * Convert from longitude, latitude, altitude (above geoid) to spherical polar coordinates
 *
 * Finishes the job started in AnitaGeomTool
 * 
 * @param lon is the longitude
 * @param lat is the latitude 
 * @param alt is the altitude above the geoid in metres
 * @param r is the radial position is metres
 * @param phi is the the azimuthal angle (increasing east/(west?) from Greenwich meridian)
 * @param theta is the elevation angle (theta = 0 points to north, theta = pi points south)
 */
void lonLatAltToSpherical(double lon, double lat, double alt, double& r, double& theta, double& phi){
  double cartesian[3];
  AnitaGeomTool::Instance()->getCartesianCoords(lat, lon, alt, cartesian);
  double x = cartesian[0];
  double y = cartesian[1];
  double z = cartesian[2];

  // AnitaGeomTool confuses and infuriates in equal parts...
  z = lat >= 0 ? TMath::Abs(z) : -TMath::Abs(z);

  r = TMath::Sqrt(x*x + y*y + z*z);
  theta = r > 0 ? TMath::ACos(z/r) : 0;
  phi = -TMath::ATan2(y, x) + 0.5*TMath::Pi();
  phi = phi >= TMath::Pi() ?  phi - TMath::TwoPi() : phi;

  // std::cout << lon <<  "\t" << lat << "\t" << alt  << std::endl;
  // std::cout << phi*TMath::RadToDeg() << "\t" << theta*TMath::RadToDeg() << "\t" << r << std::endl << std::endl;
}



/** 
 * Convert from spherical polar (r, theta, phi) to lon, lat, alt
 * 
 * @param lon is the longitude (degrees)
 * @param lat is the latitude (degrees)
 * @param alt is the altitude (meters)
 * @param r is the radial position (meters)
 * @param theta is the elevation angle (radians), theta = 0 at the north pole, increases to pi at the south pole
 * @param phi is the azimuthal angle (radians), east is +ve, west is -ve.
 */
void sphericalToLatLonAlt(double& lon, double& lat, double& alt, double r, double theta, double phi){

  double x = r*TMath::Sin(phi)*TMath::Sin(theta);
  double y = r*TMath::Cos(phi)*TMath::Sin(theta);
  double z = r*TMath::Cos(theta);
  double cartesian[3] = {x, y, z};

  auto g = AnitaGeomTool::Instance();
  g->getLatLonAltFromCartesian(cartesian, lat, lon, alt);

  // fml... 
  lat = theta*TMath::RadToDeg() <= 90 ? -lat : lat;
}





/** 
 * Get the g Gauss coefficient in the IGRF/DGRF model
 *
 * Interpolates the between the known values or extrapolates otherwise
 
 * @param t is the unixTime (seconds since 1970)
 * @param n is the degree
 * @param m is the order
 * 
 * @return the time interpolated IGRF/DGRF g coefficient
 */
double GeoMagnetic::g(UInt_t unixTime, int n, int m){
  init();
  int year = 2015;
  int index = getIndex(n, m);
  return g_vs_time[year].at(index);
}


/** 
 * Get the h Gauss coefficient in the IGRF/DGRF model
 *
 * Interpolates the between the known values or extrapolates otherwise
 
 * @param t is the unixTime (seconds since 1970)
 * @param n is the degree
 * @param m is the order
 * 
 * @return the time interpolated IGRF h coefficient
 */
double GeoMagnetic::h(UInt_t unixTime, int n, int m){
  init();  
  int year = 2015;
  int index = getIndex(n, m);  
  return h_vs_time[year].at(index);
}










/** 
 * Evaluate the Schmidt Quazi normalized asssociated Legendre polynomal at x
 *
 * Thankfully the Associated Legendre Polynomials are implemented inside ROOT
 * (actually ROOT just wraps the GSL implementation)
 * 
 * @param n is the degree valid for n>0
 * @param m is the order valid m=0, m<=n
 * @param x is the value at which to evaluate the polynomial (valid x>=-1, x<=1)
 * 
 * @return the value of the polynomial
 */
double evalSchmidtQuasiNormalisedAssociatedLegendre(int n, int m, double x){
 
  double norm = m == 0 ? 1 : TMath::Sqrt(2*getFactorial(n-m)/getFactorial(n+m));
  double P_n_m = norm*fAssocLegendre[n][m]->Eval(x);

  if(TMath::IsNaN(P_n_m)){
    std::cerr << "You got a NaN... " << norm << "\t"  << m << "\t" << n << "\t" << x << std::endl;
  }  
  return P_n_m;
}

/** 
 * @brief Workhorse function for calculating the geomagnetic potential at any point above the earth spherical polar coordinates
 * Uses the GeoMagnetic model.
 *
 * @param r is the radius from the Earth centre in metres
 * @param theta is the angle from the equator in radians, north is +ve
 * @param phi is the angle from the Greenwich maridian (lon=0) in radians
 * @param t is the time, currently just a dummy variable...
 * 
 * @return 
 */
double GeoMagnetic::getPotentialAtSpherical(UInt_t unixTime, double r, double theta, double phi){
  init();
  int year = 2015; // for now,  should be a function of time  
  double V = 0; // the potential

  // sum over the legendre polynomials normalised by the Gauss coefficients g and h
  for(int n=1;  n < numPoly; n++){
    for(int m=0;  m <= n;  m++){
      double mPhi = m*phi;
      double part = 0;

      double this_g = GeoMagnetic::g(unixTime, n, m);
      if(this_g != 0){
        part += this_g*TMath::Cos(mPhi);
      }
      double this_h = GeoMagnetic::h(year, n, m);
      if(this_h){
        part += this_h*TMath::Sin(mPhi);
      }

      if(part != 0){
        double cosTheta = TMath::Cos(theta);
        double P_n_m = evalSchmidtQuasiNormalisedAssociatedLegendre(n, m, cosTheta);
        part *= earth_radius*pow(earth_radius/r, n+1)*P_n_m;
        V += part;

        if(TMath::IsNaN(part)){
          std::cout << earth_radius << "\t" << r << "\t" << pow(earth_radius/r, n+1) << "\t" << this_g << "\t" << cos(mPhi) << "\t" <<  this_h << "\t" << sin(mPhi) << "\t" << P_n_m << std::endl;
        }
      }
    }
  }
  
  return V;
}










/** 
 * Calculates the X, Y, Z components of the geo-magnetic field
 * 
 * @param unixTime 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return a Field object containing the vector field components
 */
GeoMagnetic::Field GeoMagnetic::getFieldAtLonLatAlt(UInt_t unixTime, double lon, double lat, double alt){
  init();
  double r, theta, phi;
  lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  return getFieldAtSpherical(unixTime, r, theta, phi);
}



/** 
 * Calculates the X, Y, Z components of the geo-magnetic field
 * 
 * @param unixTime 
 * @param r 
 * @param theta 
 * @param phi 
 * 
 * @return a Field object containing the vector field components
 */
GeoMagnetic::Field GeoMagnetic::getFieldAtSpherical(UInt_t unixTime, double r, double theta, double phi){
  init();

  // There are some efficiency gains to be had if anyone wants them.
  // Each of these functions calculates the potential at r, theta, phi (and each a small, orthogonal distance away)
  // This could be done just once to save 2/6 geomagnetic potnential calculations...
  
  Field f;
  f.X = X_atSpherical(unixTime, r,  theta, phi);
  f.Y = Y_atSpherical(unixTime, r,  theta, phi);
  f.Z = Z_atSpherical(unixTime, r,  theta, phi);

  // since the direction of the X, Y, Z changes with position
  // I will put the source location in the back end of the vector
  f.r = r;
  f.theta = theta;
  f.phi = phi;
  return f;
}


/** 
 * @brief Get the geomagnetic potential at a particular time,  longitude, latitude, and altitude
 * 
 * @param unixTime 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return 
 */
double GeoMagnetic::getPotentialAtLonLatAlt(UInt_t unixTime, double lon, double lat, double alt){
  init();
  
  double r, theta, phi;
  lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  return getPotentialAtSpherical(unixTime, r, theta, phi);
}




/** 
 * Get the northwards component of the geo-magnetic field,
 * 
 * 
 * @param unixTime 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return 
 */
double GeoMagnetic::X_atLonLatAlt(UInt_t unixTime, double lon, double lat, double alt){  
  init();
  double r, phi, theta;
  lonLatAltToSpherical(lon, lat, alt, r, theta, phi);

  return X_atSpherical(unixTime, r, theta, phi);
}

/** 
 * Get the northwards component of the geo-magnetic field,
 * 
 * 
 * @param unixTime 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return 
 */


double GeoMagnetic::X_atSpherical(UInt_t unixTime, double r, double theta, double phi){
// double GeoMagnetic::X(UInt_t unixTime, double lon, double lat, double alt){  
  
  init();
//   double r, phi, theta;
//   lonLatAltToSpherical(lon, lat, alt, r, theta, phi);

  double V0 = getPotentialAtSpherical(unixTime, r, theta, phi);
  double V1 = getPotentialAtSpherical(unixTime, r, theta+dTheta, phi);
  double BX = (V1-V0)/(dTheta*r);
  return BX;
}




/** 
 * Get the east/(west?) component of the geo-magnetic field,
 * 
 * @param unixTime 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return 
 */
double GeoMagnetic::Y_atLonLatAlt(UInt_t unixTime, double lon,  double lat, double alt){
  init();
  double r, phi, theta;
  lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  return Y_atSpherical(unixTime, lon, lat, alt);
}


/** 
 * Get the east/(west?) component of the geo-magnetic field,
 * 
 * @param unixTime 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return 
 */
double GeoMagnetic::Y_atSpherical(UInt_t unixTime, double r,  double theta, double phi){
  init();
  double V0 = getPotentialAtSpherical(unixTime, r, theta, phi);
  double V1 = getPotentialAtSpherical(unixTime, r, theta, phi+dPhi);
  double BY = -(V1-V0)/(dPhi*r*TMath::Sin(theta));
  return BY;
}




/** 
 * Get the downwards facing component of the geomagnetic fielda
 * 
 * @param unixTime 
 * @param lon is the longitude, -ve is east, +ve is west (degrees)
 * @param lat is latitude, +ve is north, -ve is south
 * @param alt is altitude above geoid surface
 * 
 * @return Downwards component of geom-magnetic field
 */
double GeoMagnetic::Z_atLonLatAlt(UInt_t unixTime, double lon, double lat, double alt){  
  init();
  double r, theta, phi;
  lonLatAltToSpherical(lon, lat, alt, r, theta,  phi);
  return Z_atSpherical(unixTime, r, theta, phi);
}


/** 
 * Get the downwards facing component of the geomagnetic field
 * 
 * @param unixTime 
 * @param r 
 * @param theta 
 * @param phi 
 * @param double 
 * 
 * @return Downwards component of geomagnetic field
 */
double GeoMagnetic::Z_atSpherical(UInt_t unixTime, double r,  double theta, double phi){
  init();
  
  double V0 = getPotentialAtSpherical(unixTime, r, theta, phi);
  double V1 = getPotentialAtSpherical(unixTime, r+dr, theta, phi);
  double BZ = -(V1-V0)/dr; // negative of the gradient of the potential
  return BZ;
}


