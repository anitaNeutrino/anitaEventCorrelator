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
#include "UsefulAdu5Pat.h"
#include "TGraphAntarctica.h"

#include "TDatime.h"
#include "TCanvas.h"

// --------------------------------------------------------------------------------------------------------------------------------------
// Silly globals, kept tucked away from prying eyes
// --------------------------------------------------------------------------------------------------------------------------------------

bool doneInit = false; // Tells you whether we've read in the data and precalculated the factorials
const int numPoly = 14; // there are only 13 polynomial coeffients, but I'm going to start counting from one for simplicity
std::vector<double> factorials(2*numPoly, 0);
std::map<int, std::vector<double> > g_vs_time; // Gauss coefficients (needed to calc potential)
std::map<int, std::vector<double> > h_vs_time; // Gauss coefficients (needed to calc potential)
const double earth_radius = 6371.2e3; // earth radius in meters for the magnetic model
TF1* fAssocLegendre[numPoly][numPoly] = {{NULL}}; // Associated Legendre polynomials
bool debug = false;




/** 
 * Will print a bunch of stuff to the screen and draw pretty-ish pictures
 * 
 * @param db debugging boolian
 */
void GeoMagnetic::setDebug(bool db){
  debug = db;
}

// for differentiating the potential
const double dr = 1;
const double dTheta = 0.01*TMath::DegToRad();  
const double dPhi = 0.01*TMath::DegToRad();


// for the atmospheric model
TGraph grAtmosDensity;
TF1* fExpAtmos;
const double xMax = 0.8e4; // kg / m^{2}


// --------------------------------------------------------------------------------------------------------------------------------------
// Utility functions, for initialisation and internal calculation
// --------------------------------------------------------------------------------------------------------------------------------------

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

    const int numAtmospherePoints = 8;    
    double heights[numAtmospherePoints]   = {-611,   11019,  20063,  32162,  47350,  51413, 71802, 86000};
    double densities[numAtmospherePoints] = {1.2985, 0.3639, 0.0880, 0.0132, 0.0020, 0,     0,     0    };
    for(int i=0;  i < numAtmospherePoints; i++){
      grAtmosDensity.SetPoint(i, heights[i], densities[i]);      
    }
    grAtmosDensity.SetTitle("Atmospheric density (International Standard); Altitude above MSL (m); Densities (kg/m^{3})");
    grAtmosDensity.SetName("grAtmosDensity");    
    grAtmosDensity.Fit("expo", "Q");
    fExpAtmos = (TF1*) grAtmosDensity.FindObject("expo");
    if(!fExpAtmos){
      std::cerr << "Warning in " << __FILE__ << " unable to find exponential fit to atmosphere..." << std::endl;
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
 * Lon lat alt to cartesian TVector3
 * 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return 
 */
TVector3 lonLatAltToVector(double lon, double lat, double alt){
  double r, theta, phi;
  lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  TVector3 v;
  v.SetMagThetaPhi(r, theta, phi);
  return v;
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
 * Convert from a cartesian TVector3 to lon, lat, alt
 * 
 * @param lon 
 * @param lat 
 * @param alt 
 * @param v 
 */
void vectorToLonLatAlt(double& lon, double& lat, double& alt, const TVector3& v){
  sphericalToLatLonAlt(lon, lat, alt, v.Mag(), v.Theta(), v.Phi());
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
  init();
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
 * Plots arrows representing the B field direction to visualise the magnetic field over Antarctica
 * 
 * @param unixTime is the time
 * @param altitude is the altitude at which to calculate the B-field
 * 
 * @return the canvas on which the plot is produced
 */
TCanvas* GeoMagnetic::plotFieldAtAltitude(UInt_t unixTime, double altitude){

  auto c = new TCanvas();
  AntarcticaBackground* bg = new AntarcticaBackground();

  int nx = bg->GetNbinsX();
  int ny = bg->GetNbinsY();
  bg->Draw();

  const int arrowEvery = 20;
  for(int by=1; by <= ny; by+=arrowEvery){
    double northing = bg->GetYaxis()->GetBinLowEdge(by);
    for(int bx=1; bx <= nx; bx+=arrowEvery){
      double easting = bg->GetXaxis()->GetBinLowEdge(bx);
      double lon, lat;
      RampdemReader::EastingNorthingToLonLat(easting, northing, lon, lat);
      FieldPoint* f = new FieldPoint(0, lon, lat, 0);
      f->SetBit(kMustCleanup);
      f->SetBit(kCanDelete);
      f->Draw();
    }
  }
  return c;
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
  double BZ = (V1-V0)/dr; // negative of the gradient of the potential
  return BZ;
}




/** 
 * Handles conversion to easting/northing for use with AntarcticaBackground
 * 
 * @param opt the draw option, passed to the TArrow draw option
 */
void GeoMagnetic::FieldPoint::Draw(Option_t* opt){

  double lon, lat, alt;
  double r = fPosition.Mag();
  double theta = fPosition.Theta();
  double phi = fPosition.Phi();

  sphericalToLatLonAlt(lon, lat, alt, r, theta, phi);
  
  RampdemReader::LonLatToEastingNorthing(lon, lat, fX1, fY1);

  TVector3 position = fPosition;
  position += fDrawScaleFactor*fField;
  sphericalToLatLonAlt(lon, lat, alt, position.Mag(), position.Theta(), position.Phi());
  RampdemReader::LonLatToEastingNorthing(lon, lat, fX2, fY2);
  
  TArrow::Draw(opt);
}



GeoMagnetic::FieldPoint::FieldPoint(UInt_t unixTime, double lon, double lat, double alt) : TArrow(0, 0, 0, 0, 0.001, "|>"), fPosition(),fField(), fDrawScaleFactor(10)
{
  double r, theta, phi;
  lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  fPosition.SetMagThetaPhi(r, theta, phi);

  // each of these functions calcuates calculates V0, so you could save 2 of 6 calculations here...
  double X = X_atSpherical(unixTime, r,  theta, phi);
  double Y = Y_atSpherical(unixTime, r,  theta, phi);
  double Z = Z_atSpherical(unixTime, r,  theta, phi);

  // now convert into proper spherical polar coordinates..
  
  // X points north, but thetaHat increases to the south, so *= -1
  double B_theta = -X;
  // Y points east, if you add positive values of phi, you're going east, which is the same as spherical coordinates so the sign is the same
  double B_phi = Y;
  // Z points down, but we want the radial component to point away from the origin, so *=-1
  double B_r = -Z;

  // now I'm going to rotate into a cartesian coordinate system with the +ve z-axis running through the geographic north pole  
  double cos_theta = TMath::Cos(theta);
  double sin_theta = TMath::Sin(theta);
  double cos_phi = TMath::Cos(phi);
  double sin_phi = TMath::Sin(phi);

  // this rotates the field components pointing along RHat, thetaHat, phiHat into the Cartesian coordinate system
  double x = sin_theta*cos_phi*B_r + cos_theta*cos_phi*B_theta - sin_phi*B_phi;
  double y = sin_theta*sin_phi*B_r + cos_theta*sin_phi*B_theta + cos_phi*B_phi;
  double z = cos_theta*B_r         - sin_theta*B_theta         + 0;

  fField.SetXYZ(x, y, z);  
}



/** 
 * Does some specular reflection, make sure the "incident vector" points FROM the source to the reflection point
 * 
 * @param reflectionPointToSource 
 * @param surfaceNormal 
 * 
 * @return reflected vector
 */
TVector3 GeoMagnetic::reflection(const TVector3& sourceToReflection, const TVector3& surfaceNormal){

  // https://en.wikipedia.org/wiki/Specular_reflection#Vector_formulation
  TVector3 reflectionPointToSource = -1*sourceToReflection;
  TVector3 reflectionPointToDestination = 2*(reflectionPointToSource.Dot(surfaceNormal))*surfaceNormal - reflectionPointToSource;
  return reflectionPointToDestination;
}



/** 
 * Get a unit length TVector that points along thetaWave and phiWave 
 * 
 * @param usefulPat is ANITA's position
 * @param phiWave is the incoming azimuth direction (radians) in payload coordinates
 * @param thetaWave is the elevation angle (radians) theta=0 lies along the horizonal with -ve theta being up (the UsefulAdu5Pat convention)
 * 
 * @return TVector3 containing a unit vector pointing to thetaWave/phiWave away from ANITA
 */

TVector3 GeoMagnetic::getUnitVectorAlongThetaWavePhiWave(UsefulAdu5Pat& usefulPat, double phiWave, double thetaWave){
  init();
  
  TVector3 anitaPosition = lonLatAltToVector(usefulPat.longitude, usefulPat.latitude, usefulPat.altitude);
  
  // now I need to get a vector pointing along thetaWave and phiWave from ANITA's position
  // so let's get theta and phi wave from an arbitrary position close to the payload,
  // evaluate theta/phi expected and rotate that vector around ANITA's position
  // until it aligns with phiWave...

  // This is just due north of ANITA
  double testLon = usefulPat.longitude;
  double testLat = usefulPat.latitude - 0.1; // if ANITA could be at the north pole, this wouldn't work
  double testAlt = usefulPat.altitude;

  double testThetaWave, testPhiWave;
  usefulPat.getThetaAndPhiWave(testLon, testLat, testAlt, testThetaWave, testPhiWave);

  TVector3 testVector = lonLatAltToVector(testLon, testLat, testAlt);

  if(debug){
    std::cout << "Before rotation..." << std::endl;
    std::cout << testLon << "\t" << testLat << "\t" << testAlt << std::endl;
    std::cout << testThetaWave << "\t" << testPhiWave << std::endl;
    std::cout << testThetaWave*TMath::RadToDeg() << "\t" << testPhiWave*TMath::RadToDeg() << std::endl;
  }
  
  // payload phi increases anticlockwise (around +ve z axis)
  // the TVector3 phi increases clockwise (around +ve z axis)

  const TVector3 unitAnita = anitaPosition.Unit();
  testVector.Rotate(-testPhiWave, unitAnita); // if we were to recalculate the phiWave expected, it would now point to 0
  testVector.Rotate(phiWave, unitAnita); // if we were to recalculate the phiWave expected, it would now point to phiWave

  vectorToLonLatAlt(testLon, testLat, testAlt, testVector);
  usefulPat.getThetaAndPhiWave(testLon, testLat, testAlt, testThetaWave, testPhiWave);
    
  if(debug){
    std::cout << "After phi rotation..." << std::endl;
    std::cout << testLon << "\t" << testLat << "\t" << testAlt << std::endl;
    std::cout << testThetaWave << "\t" << testPhiWave << std::endl;
    std::cout << testThetaWave*TMath::RadToDeg() << "\t" << testPhiWave*TMath::RadToDeg() << std::endl;
  }
  
  // now need to raise/lower the point described by testVector such that thetaWave is correct
  // i.e. set the magnitude of the testVector such that thetaWave is correct
  //
  //                        
  //                      O   
  // Earth Centre (Origin) o
  //                       |\  a 
  //                     t | \ 
  //                       |  \
  //                 ANITA o---o "Test Vector"
  //                       A    T
    
  // O, A, T are the angles
  // o, a, t are the lengths. I'm trying to find the length a for a given angle A.
  // a / sin(A) = t / sin(T)
  //
  // t = anitaPosition.Mag();
  // A = (pi/2 - thetaWave);
  // O = angle between ANITA and the test vector
  // T = pi - A - O
  // so...
  // a = sin(A) * t / (pi - A - O)
  double A = TMath::PiOver2() - thetaWave;
  double O = testVector.Angle(anitaPosition); //angle between the vectors
  double T = TMath::Pi() - A - O;
  double t = anitaPosition.Mag();
  double a = TMath::Sin(A)*(t/TMath::Sin(T));

  if(debug){
    std::cout << thetaWave*TMath::RadToDeg() << "\t" << A*TMath::RadToDeg() << "\t" << O*TMath::RadToDeg() << "\t" << T*TMath::RadToDeg() << std::endl;
    std::cout << a << "\t" << t << "\t" << a - t << std::endl;
  }
  
  testVector.SetMag(a);

  if(debug){
    vectorToLonLatAlt(testLon, testLat, testAlt, testVector);
    usefulPat.getThetaAndPhiWave(testLon, testLat, testAlt, testThetaWave, testPhiWave);
    std::cout << "After setting mangitude..." << std::endl;
    std::cout << testLon << "\t" << testLat << "\t" << testAlt << std::endl;
    std::cout << testThetaWave << "\t" << testPhiWave << std::endl;
    std::cout << testThetaWave*TMath::RadToDeg() << "\t" << testPhiWave*TMath::RadToDeg() << std::endl;
  }

  TVector3 deltaVec = testVector - anitaPosition;
  return deltaVec.Unit();
}





/** 
 * Gets the expected polarisation angle for a 1e19 eV cosmic ray traversing the atmosphere
 * 
 * @param usefulPat contains ANITA's position
 * @param phiWave azimith
 * @param thetaWave elevation angle (radians) in payload coordinates (-ve theta is up, +ve theta is down)
 * 
 * @return 
 */

double GeoMagnetic::getExpectedPolarisation(UsefulAdu5Pat& usefulPat, double phiWave, double thetaWave){

  init();
  
  if(debug){
    std::cout << "Anita position = " << usefulPat.longitude << "\t" << usefulPat.latitude << "\t" << usefulPat.altitude << std::endl;
  }
  
  // use the silly UsefulAdu5Pat convention that -ve theta is down...
  // phiWave is in radians relative to ADU5 Aft Fore line

  double reflectionLon=0, reflectionLat=0, reflectionAlt=0, deltaTheta=100; // need non-zero deltaTheta when testing whether things intersectg, as theta < 0 returns instantly
  usefulPat.traceBackToContinent(phiWave, thetaWave, &reflectionLat, &reflectionLon, &reflectionAlt, &deltaTheta);

  TVector3 destination; // ANITA position if direct cosmic ray or surface position if reflected cosmic ray
  TVector3 destinationToSource; // unit vector from ANITA (if direct) or reflection position (if indirect) which points in the direction the cosmic ray signal came from
  bool directCosmicRay = TMath::Abs(deltaTheta) > 1e-20 ? true : false; // direct cosmic ray?

  TVector3 anitaPosition = lonLatAltToVector(usefulPat.longitude, usefulPat.latitude, usefulPat.altitude);

  // here we do the geometry slightly differently for the direct vs. reflected case
  if(directCosmicRay){
    if(debug){
      std::cout << "I'm a direct Cosmic Ray! " << deltaTheta << "\t" <<  reflectionLon << "\t" << reflectionLat << "\t"  << reflectionAlt  << std::endl;
    }
    destination = anitaPosition;
    destinationToSource = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave, thetaWave);
  }
  else{ // reflected cosmic ray
    if(debug){
      std::cout << "I'm a reflected Cosmic Ray! " << deltaTheta << "\t" <<  reflectionLon << "\t" << reflectionLat << "\t"  << reflectionAlt  << std::endl;
    }
    destination = lonLatAltToVector(reflectionLon, reflectionLat, reflectionAlt); // i.e. the reflection point

    TVector3 reflectionToAnita = anitaPosition - destination; // from reflection to anita...

    // here I find the normal to  the geoid surface by getting the vector difference between
    // a point 1 m above the reflection and the reflection
    TVector3 surfaceNormal = (lonLatAltToVector(reflectionLon, reflectionLat, reflectionAlt + 1) - destination).Unit();

    // Reflect the incoming vector...
    TVector3 incomingVector = reflection(reflectionToAnita, surfaceNormal);

    // but I want the vector pointing from the reflection out towards where the income came from
    destinationToSource = -incomingVector.Unit();
  }

  // Here we get the position the cosmic way was on top of the atmosphere...
  TVector3 topOfAtmosphere = getInitialPosition(destination, destinationToSource);

  // ...And the direction it travelled through the atmosphere...
  const TVector3 cosmicRayDirection = -destinationToSource.Unit();

  // We integrate along that vector with our atmospheric model until we get to xMax for the shower
  TVector3 xMaxPosition = getXMaxPosition(topOfAtmosphere, cosmicRayDirection);

  // and calculate the geo-magnetic field field at the shower maximum
  double xMaxLon=0, xMaxLat=0, xMaxAlt=0;
  vectorToLonLatAlt(xMaxLon, xMaxLat, xMaxAlt, xMaxPosition);
  FieldPoint fp(usefulPat.realTime, xMaxLon, xMaxLat, xMaxAlt);
  
  // This is our electric field vector!
  // If we care about getting the magnitude correct in addition to the orientation, there are some missing factors
  // It should be:
  // B_vec x S_vec = (1/mu0)B^{2} E_vec
  // but for now this will do.
  TVector3 EVec = fp.field().Cross(cosmicRayDirection).Unit();
  
  if(!directCosmicRay){
    std::cerr << "Need to implement E-Field reflection with Fresnel coefficients and other bells and whistles." << std::endl;
  }

  // Here I find the VPol and HPol antenna axes.
  // And I'm going to  pretend that one of ANITA's antennas points exactly at phiWave

  // Since the antennas points down at -10 degrees, the VPol axis is 80 degrees above the horizontal plane
  TVector3 vPolAxis = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave, -80*TMath::DegToRad());
  // The VPol feed is up... (if) the HPol feed is to the right (looking down the boresight) then it points anticlockwise around the payload
  // phi increases anti-clockwise in payload coordinates, therefore
  TVector3 hPolAxis = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave + TMath::PiOver2(), -10*TMath::DegToRad());
  
  // Dot the electric field with the antenna polarisation vectors...
  double vPolComponent = EVec.Dot(vPolAxis);
  double hPolComponent = EVec.Dot(hPolAxis);
  
  // et voila
  double polarisationAngle = TMath::ATan2(vPolComponent, hPolComponent);
  
  return polarisationAngle;
}












TCanvas* GeoMagnetic::plotAtmosphere(){
  init();
  TCanvas* c1 = new TCanvas();
  grAtmosDensity.Draw("al");
  return c1;
}



double GeoMagnetic::getAtmosphericDensity(double altitude){
  init();
  // TODO, convert altitude (above geoid) to height above MSL...
  return fExpAtmos->Eval(altitude);
}





TVector3 GeoMagnetic::getInitialPosition(const TVector3& destination, const TVector3& destinationToSource){
  init();

  if(destinationToSource.Mag() != 1){
    std::cerr  << "Warning in " << __PRETTY_FUNCTION__ << ", was expecting a unit vector, didn't get one. "
               << "This calculation might be wonky... " << std::endl;
  }
  
  
  // moves out from the destination in 1km steps until the altitude is 80km
  const double desiredAlt = 80e3; // 100 km
  TVector3 initialPosition = destination;
  double initialLon, initialLat, initialAlt;
  vectorToLonLatAlt(initialLon, initialLat, initialAlt, initialPosition);
  if(debug){
    std::cout  <<  "Initial alt = " << initialAlt << ", (desiredAlt = "  << desiredAlt << ")" << std::endl;
  }
  while(initialAlt < desiredAlt){
    initialPosition += 1e3*destinationToSource;
    vectorToLonLatAlt(initialLon, initialLat, initialAlt, initialPosition);
    // std::cout << initialLon << "\t" << initialLat << "\t" << initialAlt << std::endl;
  }
  if(debug){
    std::cout << "Got initial position... " << initialLon << "\t" << initialLat << "\t" << initialAlt << std::endl;    
  }
  
  return initialPosition;
}


TVector3 GeoMagnetic::getXMaxPosition(const TVector3& initialPosition, const TVector3& cosmicRayDirection){
  init();

  TVector3 currentPosition = initialPosition;
  double currentAtmosphereTraversed = 0; // kg / m ^{2}
  double dx = cosmicRayDirection.Mag();

  int numSteps = 0;

  TGraph* grAltPath = debug ? new TGraph() : NULL;

  double lastAlt = 0;
  double initialAlt = 0;
  while(currentAtmosphereTraversed < xMax){
    
    currentPosition += cosmicRayDirection;
    double currentLon, currentLat, currentAlt;
    vectorToLonLatAlt(currentLon, currentLat, currentAlt, currentPosition);

    if(initialAlt==0){
      initialAlt= currentAlt;
    }


    if(debug && currentAtmosphereTraversed == 0){
      std::cout << "Finding xMax position... Initial position... " << currentLon << "\t" << currentLat << "\t" << currentAlt << "\t" << currentAtmosphereTraversed << std::endl;
    }
    
    currentAtmosphereTraversed += getAtmosphericDensity(currentAlt)*dx;

    if(debug && currentAtmosphereTraversed >= xMax){
      std::cout << "Found xMax position... " << currentLon << "\t" << currentLat << "\t" << currentAlt << "\t" << currentAtmosphereTraversed << std::endl;
    }

    // then we're ascending without reaching xMax, which is bad
    if(currentAlt > lastAlt && lastAlt > initialAlt){
      std::cerr << "Warning in "  <<__PRETTY_FUNCTION__ << ", unable to find xMax, terminating loop." << std::endl;
      std::cerr << "Current position = " << currentLon << "\t" << currentLat << "\t" << currentAlt << "\t" << currentAtmosphereTraversed << std::endl;
      break;
    }

    if(debug){
      grAltPath->SetPoint(numSteps, dx*numSteps, currentAlt);
      numSteps++;
    }
    lastAlt = currentAlt;
  }

  if(debug){
    TCanvas* c1 = new TCanvas();
    grAltPath->SetTitle("Cosmic Ray Altitude vs. distance traversed; Distance through atmosphere (m); Altitude (m)");
    grAltPath->Draw("al");
  }
  
  return currentPosition;
}
