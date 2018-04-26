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
#include "TLegend.h"

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

// for "differentiating" the potential (really it's a difference with small delta values)
// I suppose I could differentiate the expression by hand and evaluate that... but that's work
const double dr = 1;
const double dTheta = 0.01*TMath::DegToRad();  
const double dPhi = 0.01*TMath::DegToRad();


// Globals for the atmospheric model, currently just an exponential fit to some numbers
// from a table in the Wikipedia
TGraph grAtmosDensity;
TF1* fExpAtmos;

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
void prepareGeoMagnetics(){

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
 * @param phi is the the azimuthal angle (increasing east from Greenwich meridian)
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

  AnitaGeomTool* g = AnitaGeomTool::Instance();
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










// --------------------------------------------------------------------------------------------------------------------------------------
// Namespace functions, for public consumption
// --------------------------------------------------------------------------------------------------------------------------------------




/** 
 * @brief 3D geometry is hard, set this true to get more info.
 *
 * Setting this flag will print a bunch of stuff to the screen and maybe draw pretty-ish plots as functions are called
 * 
 * @param db debugging boolian
 */
void GeoMagnetic::setDebug(bool db){
  debug = db;
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
  prepareGeoMagnetics();
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
  prepareGeoMagnetics();  
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
  prepareGeoMagnetics();
  int year = 2015; // for now, should be a function of time
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
  prepareGeoMagnetics();
  
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
  prepareGeoMagnetics();
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
  prepareGeoMagnetics();
  double V0 = getPotentialAtSpherical(unixTime, r, theta, phi);
  double V1 = getPotentialAtSpherical(unixTime, r, theta+dTheta, phi);
  double BX = (V1-V0)/(dTheta*r);
  return BX;
}




/** 
 * Get the eastwards component of the geo-magnetic field,
 * 
 * @param unixTime 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return 
 */
double GeoMagnetic::Y_atLonLatAlt(UInt_t unixTime, double lon,  double lat, double alt){
  prepareGeoMagnetics();
  double r, phi, theta;
  lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  return Y_atSpherical(unixTime, lon, lat, alt);
}


/** 
 * Get the eastwards component of the geo-magnetic field,
 * 
 * @param unixTime 
 * @param lon 
 * @param lat 
 * @param alt 
 * 
 * @return 
 */
double GeoMagnetic::Y_atSpherical(UInt_t unixTime, double r,  double theta, double phi){
  prepareGeoMagnetics();
  double V0 = getPotentialAtSpherical(unixTime, r, theta, phi);
  double V1 = getPotentialAtSpherical(unixTime, r, theta, phi+dPhi);
  double BY = -(V1-V0)/(dPhi*r*TMath::Sin(theta));
  return BY;
}




/** 
 * Get the downwards facing component of the geomagnetic fielda
 * 
 * @param unixTime 
 * @param lon is the longitude, +ve is east, -ve is west (degrees)
 * @param lat is latitude, +ve is north, -ve is south
 * @param alt is altitude above geoid surface
 * 
 * @return Downwards component of geom-magnetic field
 */
double GeoMagnetic::Z_atLonLatAlt(UInt_t unixTime, double lon, double lat, double alt){  
  prepareGeoMagnetics();
  double r, theta, phi;
  lonLatAltToSpherical(lon, lat, alt, r, theta,  phi);
  return Z_atSpherical(unixTime, r, theta, phi);
}



/** 
 * Plots arrows representing the B field direction to visualise the magnetic field over Antarctica
 * 
 * The graphical objects have the kCanDelete bit set, so if you delete the canvas they are destroyed
 * 
 * @param unixTime is the time
 * @param altitude is the altitude at which to calculate the B-field
 * 
 * @return the canvas on which the plot is produced
 */
TCanvas* GeoMagnetic::plotFieldAtAltitude(UInt_t unixTime, double altitude){

  TCanvas* c = new TCanvas();
  AntarcticaBackground* bg = new AntarcticaBackground();

  int nx = bg->GetNbinsX();
  int ny = bg->GetNbinsY();
  bg->Draw();
  bg->SetBit(kCanDelete);
  
  const int arrowEvery = 20;
  for(int by=1; by <= ny; by+=arrowEvery){
    double northing = bg->GetYaxis()->GetBinLowEdge(by);
    for(int bx=1; bx <= nx; bx+=arrowEvery){
      double easting = bg->GetXaxis()->GetBinLowEdge(bx);
      double lon, lat;
      RampdemReader::EastingNorthingToLonLat(easting, northing, lon, lat);
      FieldPoint* f = new FieldPoint(unixTime, lon, lat, altitude);
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
  prepareGeoMagnetics();
  
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






/** 
 * Constructor for the field point taking a TVector containing the cartesian position at which you want to evaluate the field
 * 
 * @param unixTime is the time at which you wish to evaluate the field
 * @param position is the cartesian position at which you wish to evaluate the field
 */
GeoMagnetic::FieldPoint::FieldPoint(UInt_t unixTime, const TVector3& position) : TArrow(0, 0, 0, 0, 0.001, "|>"), fPosition(),fField(), fDrawScaleFactor(10)
{
  fPosition = position;
  fUnixTime = unixTime;
  calculateFieldAtPosition();  
}



/** 
 * Human friendly constructor for the field point taking lontitude, latitude and altitude at which you want to evaluate the field
 * 
 * @param unixTime is the time at which you wish to evaluate the field
 * @param lon is the longitude
 * @param lat is the latitude
 * @param alt is the altitude above the geoid surface
 */
GeoMagnetic::FieldPoint::FieldPoint(UInt_t unixTime, double lon, double lat, double alt) : TArrow(0, 0, 0, 0, 0.001, "|>"), fPosition(),fField(), fDrawScaleFactor(10)
{
  double r, theta, phi;
  lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  fPosition.SetMagThetaPhi(r, theta, phi);
  fUnixTime = unixTime;
  calculateFieldAtPosition();
}




/** 
 * Workhorse function to calcuate the magnetic field from the potential at fPosition
 * Everything difficult regarding coordinate transformations with the field is in here
 */
void GeoMagnetic::FieldPoint::calculateFieldAtPosition(){
  // each of these functions calcuates calculates V0, so you could save 2 of 6 calculations here...
  double r = fPosition.Mag();
  double theta = fPosition.Theta();
  double phi = fPosition.Phi();
  
  double X = X_atSpherical(fUnixTime, r,  theta, phi); // north
  double Y = Y_atSpherical(fUnixTime, r,  theta, phi); // east
  double Z = Z_atSpherical(fUnixTime, r,  theta, phi); // down

  // now convert into proper spherical polar coordinates..
  
  // X points north, but theta_hat increases to the south, so *= -1
  double B_theta = -X;
  // Y points east, if you add positive values of phi, you're going east, which is the same as spherical coordinates so the sign is the same
  double B_phi = Y;
  // Z points down, but we want the radial component to point away from the origin, so *=-1
  double B_r = -Z;

  // Red for radially upwards pointing B-field
  // Blue for radially downwards pointing B-field
  SetLineColor(B_r > 0 ? kRed : kBlue);

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
 * Does some specular reflection, make sure the incident vector points FROM the source TO the reflection point
 * 
 * @param reflectionPointToSource incoming vector
 * @param surfaceNormal defines the plane of reflection
 * 
 * @return reflected vector
 */
TVector3 GeoMagnetic::specularReflection(const TVector3& incidentPoyntingVector, const TVector3& surfaceNormal){

  // https://en.wikipedia.org/wiki/Specular_reflection#Vector_formulation
  TVector3 reflectionPointToSource = -1*incidentPoyntingVector;
  TVector3 reflectionPointToDestination = 2*(reflectionPointToSource.Dot(surfaceNormal))*surfaceNormal - reflectionPointToSource;
  return reflectionPointToDestination;
}




/** 
 * Draw the fresnel electric field coefficients to test they are correctly implemented.
 * 
 * The graphical objects have the kCanDelete bit set, so if you delete the canvas they are destroyed
 * 
 * @return the canvas upon which the plots were drawn
 */
TCanvas* GeoMagnetic::plotFresnelReflection(){

  const int nSteps = 90; //1000;
  double d_theta_i = TMath::PiOver2()/nSteps;

  const TVector3 surfaceNormal(0, 0, 1); // reflection in the x-y plane, == z_hat
  const TVector3 x_hat(1, 0, 0);
  const TVector3 y_hat(0, 1, 0);
      
  const double pol_angle = 45*TMath::DegToRad();

  const double E_mag = 1;

  TGraph* grThetas = new TGraph();
  TGraph* gr_r_p = new TGraph();
  TGraph* gr_r_s = new TGraph();

  for(int i=0; i <= nSteps; i++){
    double theta_i = d_theta_i * i;

    TVector3 incidentRay = -surfaceNormal; // *-1 to make it downgoing
    incidentRay.RotateY(theta_i); // this should rotate the vector into the x-z plane, the plane of incidence

    TVector3 E = E_mag*y_hat; // y_hat should be perpendicular to the incident ray in the x-z plane
    E.Rotate(pol_angle, incidentRay);
   
    if(debug){
      std::cout << "theta_i = " << theta_i*TMath::RadToDeg() << " degrees, (incidentRay . E) = " << incidentRay.Dot(E) << std::endl;
      std::cout << "Incident ray = (" << incidentRay.X() << ", " << incidentRay.Y() << ", " << incidentRay.Z() << ")" << std::endl;
      std::cout << "E = (" << E.X() << ", " << E.Y() << ", " << E.Z() << ")" << std::endl;
    }

    double E_s_initial = E.Y();
    double E_p_initial = TMath::Sqrt(E.X()*E.X() + E.Z()*E.Z());

    TVector3 reflectedRay = fresnelReflection(incidentRay, surfaceNormal, E);

    if(debug){
      TVector3 npi = incidentRay.Cross(reflectedRay).Unit();
      std::cout << "normal to the plane of incidence = (" << npi.X() << ", " << npi.Y() << ", " << npi.Z() << ")" << std::endl;
    }

    double theta_r = surfaceNormal.Angle(reflectedRay);
    double theta_i_2 = surfaceNormal.Angle(-incidentRay);

    grThetas->SetPoint(grThetas->GetN(), theta_i_2*TMath::RadToDeg(), theta_r*TMath::RadToDeg());

    double E_s_reflected = E.Y();
    double E_p_reflected = TMath::Sqrt(E.X()*E.X() + E.Z()*E.Z());

    double r_p = TMath::Abs(E_p_initial) > 0 ? E_p_reflected/E_p_initial : 0; // might be a better value to choose here?
    double r_s = TMath::Abs(E_s_initial) > 0 ? E_s_reflected/E_s_initial : 0; // might be a better value to choose here?

    gr_r_p->SetPoint(gr_r_p->GetN(), theta_i_2*TMath::RadToDeg(), r_p);
    gr_r_s->SetPoint(gr_r_s->GetN(), theta_i_2*TMath::RadToDeg(), r_s);
  }    
    
  TCanvas* c1 = new TCanvas();
  c1->Divide(2);
  c1->cd(1);
  grThetas->SetTitle("Law of reflection; Angle of incidence, #theta_{i} (Degrees); Angle of reflection #theta_{r} (Degrees)");
  grThetas->Draw();
  grThetas->SetBit(kCanDelete);

  c1->cd(2);
  
  TLegend* l1b = new TLegend(0.15, 0.55, 0.6, 0.85);
  gr_r_p->SetTitle("E-field amplitude reflection coefficients; Incident angle, #theta_{i} (Degrees); Reflection coefficient");
  gr_r_p->SetLineColor(kRed);
  gr_r_p->Draw("al");
  gr_r_p->SetBit(kCanDelete);  
  double r_min = -1;
  double r_max = 1;
  gr_r_p->SetMaximum(r_max);
  gr_r_p->SetMinimum(r_min);
  gr_r_p->GetXaxis()->SetRangeUser(0, 90);
  gr_r_s->Draw("lsame");
  gr_r_s->SetBit(kCanDelete);    
  TGraph* grBrewster = new TGraph();
  double brewster_angle = TMath::RadToDeg()*TMath::ATan2(n_ice, n_air);
  grBrewster->SetPoint(0, brewster_angle, r_min);
  grBrewster->SetPoint(1, brewster_angle, r_max);
  grBrewster->SetLineStyle(2);
  grBrewster->SetLineColor(kMagenta);
  grBrewster->Draw("lsame");
  grBrewster->SetBit(kCanDelete);
  
  l1b->AddEntry(gr_r_s, "r_s (E-field perpendicular to plane of incidence)", "l");
  l1b->AddEntry(gr_r_p, "r_p (H-field perpendicular to plane of incidence)", "l");
  l1b->AddEntry(grBrewster, "Brewster angle", "l");  
  l1b->SetBit(kCanDelete);
  
  l1b->Draw();
  
  c1->cd();

  return c1;
}



/** 
 * Applies a reflection to the fresnel coefficients
 * 
 * @param incidentPoyntingVector is a vector pointing to the reflection surface 
 * @param surfaceNormal defines the plane of reflection
 * @param electricFieldVec is the electric field vector of the incident ray
 * @param n1 is the refractive index of the medium through which the incident and reflected rays are traveling (default is n_air)
 * @param n2 is the refractive index of the medium upon which the incident ray is reflected (default is n_ice)
 * 
 * @return the reflected poynting vector
 */
TVector3 GeoMagnetic::fresnelReflection(const TVector3& incidentPoyntingVector, const TVector3& surfaceNormal, TVector3& electricFieldVec, double n1, double n2){  
  
  double theta_i = surfaceNormal.Angle(-incidentPoyntingVector); // factor of -1 since the incident ray is downgoing
  double cos_theta_i = TMath::Cos(theta_i);
  double sin_theta_i = TMath::Sin(theta_i);

  // Snell's law + trig identities
  double cos_theta_t = TMath::Sqrt(1 - (n1*n1*sin_theta_i*sin_theta_i/(n2*n2)));

  const TVector3 reflectedPoyntingVector = specularReflection(incidentPoyntingVector, surfaceNormal);  

  // from the Wikipedia...
  // In this treatment, the coefficient r is the ratio of the reflected wave's complex electric field amplitude to that of the incident wave.
  // The light is split into s and p polarizations
  // s => E field is perpendicular to plane of incidence
  // p => E field is parallel to the plane of incidence
  // For s polarization, a POSITIVE r or t means that the ELECTRIC fields of the incoming and reflected (or transmitted) waves are PARALLEL, while negative means anti-parallel.
  // For p polarization, a POSITIVE r or t means that the MAGNETIC fields of the incoming and reflected (or transmitted) waves are PARALLEL, while negative means anti-parallel.

  // +ve r_s means E field (perpendicular to plane of incidence) is parallel on the way out (-ve means anti-parallel)
  // +ve r_p means H field (perpendicular to plane of incidence) is parallel on the way out (-ve means anti-parallel)
  // in the case of normal incidence, if the poynting vector is reversed but the H field is parallel, the E field must be anti-parallel...
  // therefore there should be a factor of -1 applied to the reflected electric field vector in front of r_p.
  
  TVector3 normalToPlaneOfIncidence = incidentPoyntingVector.Cross(reflectedPoyntingVector).Unit();

  if(normalToPlaneOfIncidence.Mag()==0){
    // then the incoming and outgoing vectors are anti-parallel and we are free to chose a plane of incidence.
    normalToPlaneOfIncidence.SetXYZ(0, 1, 0);
  }
  
  double r_s = (n1*cos_theta_i - n2*cos_theta_t)/(n1*cos_theta_i + n2*cos_theta_t);
  double r_p = (n2*cos_theta_i - n1*cos_theta_t)/(n2*cos_theta_i + n1*cos_theta_t);
  
  const TVector3 z_hat_local = surfaceNormal.Unit();
  const TVector3 y_hat_local = normalToPlaneOfIncidence; // therefore y_hat is definitionally perpendicular to the plane of incidence
  const TVector3 x_hat_local = y_hat_local.Cross(z_hat_local); 

  double E_s = electricFieldVec.Dot(y_hat_local);
  double E_p_x = electricFieldVec.Dot(x_hat_local);
  double E_p_z = electricFieldVec.Dot(z_hat_local);

  double reflected_E_s = r_s*E_s; // +ve r_s means E field is parallel
  double reflected_E_p_x = -r_p*E_p_x; // +ve r_p means H field is parallel, which means E-field is flipped
  double reflected_E_p_z = -r_p*E_p_z; // +ve r_p means H field is parallel, which means E-field is flipped

  TVector3 reflectedElectricFieldVec = reflected_E_p_x*x_hat_local + reflected_E_s*y_hat_local + reflected_E_p_z*z_hat_local;

  if(debug){
    std::cout << "incident angle = " << theta_i*TMath::RadToDeg() << " degrees, transmission angle = " << TMath::ACos(cos_theta_t)*TMath::RadToDeg() << " degrees" << std::endl;
    std::cout << "cos(theta_i) = " << cos_theta_i << ", cos(theta_t) = " << cos_theta_t << std::endl;
    std::cout << "r_s(theta_i) = " << r_s << ", r_p(theta_i) = " << r_p << std::endl;
    
    if(TMath::Abs((x_hat_local.Cross(y_hat_local) - z_hat_local).Mag()) > 1e-19){
      std::cout << "Checking that my unit vectors make sense..." << std::endl;
      std::cout << "|x| = " << x_hat_local.Mag() << ", |y| = " << y_hat_local.Mag() << ", |z| = " << z_hat_local.Mag() << std::endl;
      std::cout << "x.y = " << x_hat_local.Dot(y_hat_local) << ", x.z = " << x_hat_local.Dot(z_hat_local) << ", x.z = " << y_hat_local.Dot(z_hat_local) << std::endl;
      std::cout << "|(x.Cross(y)) - z| = " << (x_hat_local.Cross(y_hat_local) - z_hat_local).Mag() << std::endl;
    }
    std::cout << "Original electric field  = (" << electricFieldVec.X() << ", " << electricFieldVec.Y() << ", " << electricFieldVec.Z() << ")"  << std::endl;
    std::cout << "Reflected electric field = (" << reflectedElectricFieldVec.X() << ", " << reflectedElectricFieldVec.Y() << ", " << reflectedElectricFieldVec.Z() << ")"  << std::endl;
    std::cout << std::endl;
  }

  // modify the input electric field vector with the reflected value
  electricFieldVec = reflectedElectricFieldVec;

  // return the reflected poynting vector...
  return reflectedPoyntingVector;
}



/** 
 * Get a unit length TVector that points along thetaWave and phiWave 
 * 
 * @param usefulPat is ANITA's position
 * @param phiWave is the azimuth direction (radians) in payload coordinates
 * @param thetaWave is the elevation angle (radians) theta=0 lies along the horizonal with -ve theta being up (the UsefulAdu5Pat convention)
 * 
 * @return TVector3 containing a unit vector pointing to thetaWave/phiWave away from ANITA
 */

TVector3 GeoMagnetic::getUnitVectorAlongThetaWavePhiWave(UsefulAdu5Pat& usefulPat, double phiWave, double thetaWave){
  prepareGeoMagnetics();
  
  TVector3 anitaPosition = lonLatAltToVector(usefulPat.longitude, usefulPat.latitude, usefulPat.altitude);
  
  // now I need to get a vector pointing along thetaWave and phiWave from ANITA's position
  // so let's get theta and phi wave from an arbitrary position close to the payload,
  // evaluate theta/phi expected and rotate that vector around ANITA's position
  // until it aligns with phiWave...

  // This is just due north of ANITA
  double testLon = usefulPat.longitude;
  double testLat = usefulPat.latitude - 0.1; // if ANITA could be at the north pole, this might not work
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

double GeoMagnetic::getExpectedPolarisation(UsefulAdu5Pat& usefulPat, double phiWave, double thetaWave, 
					    double xMax){

  prepareGeoMagnetics();
  
  if(debug){
    std::cout << "Anita position = " << usefulPat.longitude << "\t" << usefulPat.latitude << "\t" << usefulPat.altitude << std::endl;
  }
  
  // use the silly UsefulAdu5Pat convention that -ve theta is down...
  // phiWave is in radians relative to ADU5 Aft Fore line

  double reflectionLon=0, reflectionLat=0, reflectionAlt=0, deltaTheta=100; // need non-zero deltaTheta when testing whether things intersectg, as theta < 0 returns instantly
  
  //histGround==1 or 2 means it hits ground, 0 means it doesn't
  int hitsGround = usefulPat.traceBackToContinent(phiWave, thetaWave, &reflectionLon, &reflectionLat, &reflectionAlt, &deltaTheta);

  TVector3 destination; // ANITA position if direct cosmic ray or surface position if reflected cosmic ray
  TVector3 destinationToSource; // unit vector from ANITA (if direct) or reflection position (if indirect) which points in the direction the cosmic ray signal came from
  bool directCosmicRay = hitsGround == 0 ? true : false; // direct cosmic ray?

  TVector3 anitaPosition = lonLatAltToVector(usefulPat.longitude, usefulPat.latitude, usefulPat.altitude);

  // Only used in the reflected case...
  TVector3 surfaceNormal;
  TVector3 reflectionToAnita;
  
  // here we do the geometry slightly differently for the direct vs. reflected case
  if(directCosmicRay){
    if(debug){
      std::cout << "I'm a direct Cosmic Ray! (" << hitsGround << ") " << deltaTheta << "\t" <<  reflectionLon << "\t" << reflectionLat << "\t"  << reflectionAlt  << std::endl;
    }
    destination = anitaPosition;
    destinationToSource = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave, thetaWave);
  }
  else{ // reflected cosmic ray
    if(debug){
      std::cout << "I'm a reflected Cosmic Ray! (" << hitsGround << ") " << deltaTheta << "\t" <<  reflectionLon << "\t" << reflectionLat << "\t"  << reflectionAlt  << std::endl;
    }
    destination = lonLatAltToVector(reflectionLon, reflectionLat, reflectionAlt); // i.e. the reflection point

    reflectionToAnita = anitaPosition - destination; // from reflection to anita...

    // here I find the normal to  the geoid surface by getting the vector difference between
    // a point 1 m above the reflection and the reflection
    surfaceNormal = (lonLatAltToVector(reflectionLon, reflectionLat, reflectionAlt + 1) - destination).Unit();

    // Reflect the incoming vector...
    TVector3 incomingVector = specularReflection(reflectionToAnita, surfaceNormal);

    // but I want the vector pointing from the reflection out towards where the incoming signal came from
    destinationToSource = -incomingVector.Unit();
  }

  // Here we get the position the cosmic ray was on top of the atmosphere...
  TVector3 cosmicRayAtmosphericEntry = getInitialPosition(destination, destinationToSource);

  // ...And the direction it travelled through the atmosphere...
  const TVector3 cosmicRayDirection = -destinationToSource.Unit();

  // We integrate along that vector with our atmospheric model until we get to xMax, where the cosmic ray would interact, on average
  TVector3 xMaxPosition = getXMaxPosition(cosmicRayAtmosphericEntry, cosmicRayDirection, xMax);

  // Calculate the geo-magnetic field field at x-max
  FieldPoint fp(usefulPat.realTime, xMaxPosition);
  if (debug){
    std::cout << "FieldPoint position:" << fp.posX() << "," << fp.posY() << ","<< fp.posZ() << std::endl; 
    std::cout << "FieldPoint vector:" << fp.componentX() << "," << fp.componentY() << ","<< fp.componentZ() << std::endl; 
  }
  // This is our electric field vector!
  // If I cared about getting the magnitude correct in addition to the orientation, there are some missing factors
  // It should be: B_vec x S_vec = (1/mu0)B^{2} E_vec
  // But there's no radio cherenkov/geomagnetic shower model or anything so for now just the cross product will do.
  TVector3 EVec = fp.field().Cross(cosmicRayDirection).Unit();
  if (debug) {
    std::cout << "EVec: (" << EVec.X() << "," << EVec.Y() << "," << EVec.Z() << ")" << std::endl;
  }
  
  if(!directCosmicRay){
    // Modifies EVec by applying the Fresnel coefficients during the reflection
    TVector3 reflectionToAnita2 = fresnelReflection(cosmicRayDirection, surfaceNormal, EVec);

    if(debug){
      double shouldBeZero = reflectionToAnita.Angle(reflectionToAnita2);
      std::cout << "Doing the reflection on the way back... the angle between reflectionToAnita and reflectionToAnita2 = " << shouldBeZero << std::endl;
    }
  }

  // Here I find the VPol and HPol antenna axes.
  // And I'm going to  pretend that one of ANITA's antennas points exactly at phiWave
  // getting an off axis response will be more complicated

  // Since the antennas points down at -10 degrees, the VPol axis is 80 degrees above the horizontal plane
  TVector3 vPolAxis = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave, 80*TMath::DegToRad());
  // The VPol feed is up... (if) the HPol feed is to the right (looking down the boresight) then it points anticlockwise around the payload
  // phi increases anti-clockwise in payload coordinates, therefore
  TVector3 hPolAxis = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave + TMath::PiOver2(), 0);//-10*TMath::DegToRad());

  if(debug){
    TVector3 antennaAxis = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave, -10*TMath::DegToRad());
    std::cout << "The dot products of any pair of the antenna axis, hPolAxis, vPolAxis should be zero..." << std::endl;
    std::cout << "antennaAxis dot hPolAxis  = " << antennaAxis.Dot(hPolAxis) << std::endl;
    std::cout << "antennaAxis dot vPolAxis  = " << antennaAxis.Dot(vPolAxis) << std::endl;
    std::cout << "hPolAxis dot vPolAxis  = " << hPolAxis.Dot(vPolAxis) << std::endl;
  }
  
  // Dot the electric field with the antenna polarisation vectors...
  double vPolComponent = EVec.Dot(vPolAxis);
  double hPolComponent = EVec.Dot(hPolAxis);
  if (debug) {
    std::cout << "vPolComponent: " << vPolComponent << std::endl;
    std::cout << "hPolComponent: " << hPolComponent << std::endl;
  }

  // et voila
  double polarisationAngle = TMath::ATan(vPolComponent/hPolComponent);
  
  return polarisationAngle;
}



/** 
 * Gets the expected polarisation angle for a 1e19 eV cosmic ray exiting the ice
 * 
 * @param usefulPat contains ANITA's position
 * @param phiWave azimith
 * @param thetaWave elevation angle (radians) in payload coordinates (-ve theta is up, +ve theta is down)
 * @param pathLength distance along the upgoing vector that the shower begins, since it isn't well defined 
 *
 * @return 
 */

double GeoMagnetic::getExpectedPolarisationUpgoing(UsefulAdu5Pat& usefulPat, double phiWave, double thetaWave,
						   double pathLength){
  if (debug) {
    std::cout << "----------------GeoMagnetic::getExpectedPolarisationUpgoing(): Begin." << std::endl;
  }
  prepareGeoMagnetics();
  

  if(debug){
    std::cout << "Anita position = " << usefulPat.longitude << "\t" << usefulPat.latitude << "\t" << usefulPat.altitude << std::endl;
  }
  
  // use the silly UsefulAdu5Pat convention that -ve theta is down...
  // phiWave is in radians relative to ADU5 Aft Fore line

  double reflectionLon=0, reflectionLat=0, reflectionAlt=0, deltaTheta=100; // need non-zero deltaTheta when testing whether things intersectg, as theta < 0 returns instantly
  int hitsGround = usefulPat.traceBackToContinent(phiWave, thetaWave, &reflectionLon, &reflectionLat, &reflectionAlt, &deltaTheta);

  TVector3 destination; // ANITA position if direct cosmic ray or surface position if reflected cosmic ray
  TVector3 destinationToSource; // unit vector from ANITA (if direct) or reflection position (if indirect) which points in the direction the cosmic ray signal came from
  bool directCosmicRay = hitsGround == 0 ? true : false; // direct cosmic ray?

  TVector3 anitaPosition = lonLatAltToVector(usefulPat.longitude, usefulPat.latitude, usefulPat.altitude);

  // here we do the geometry slightly differently for the direct vs. reflected case
  if(directCosmicRay){
    if (debug) {
      std::cout << "I'm pointed above the horizon!  I can't be an upgoing event!  Exiting..." << std::endl;
      std::cout << "deltaTheta=" << deltaTheta << " hitsGround=" << hitsGround << std::endl;
    }
    return -9999;
  }
  else{ // reflected cosmic ray
    if(debug){
      std::cout << "I'm pointed at the continent! " << deltaTheta << "\t" <<  reflectionLon << "\t" << reflectionLat << "\t"  << reflectionAlt  << std::endl;
    }
    destination = lonLatAltToVector(reflectionLon, reflectionLat, reflectionAlt); // i.e. the source on the ice
    destinationToSource = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave, thetaWave); // from ANITA to ice
  }

  const TVector3 surfaceNormal = (lonLatAltToVector(reflectionLon, reflectionLat, reflectionAlt + 1) - destination).Unit();
  
  //the cosmic ray is traveling in the opposite direction as the direction from ANITA to the ice
  const TVector3 cosmicRayDirection = -destinationToSource.Unit();

  double emergenceAngle = surfaceNormal.Angle(cosmicRayDirection);
  if (debug) {
    std::cout << "surfaceNormal: (" << surfaceNormal.X() << "," << surfaceNormal.Y() << "," << surfaceNormal.Z() << ")" << std::endl;
    std::cout << "cosmicRayDirection: (" << cosmicRayDirection.X() << "," << cosmicRayDirection.Y() << "," << cosmicRayDirection.Z() << ")" << std::endl;
    std::cout << "Emergence Angle: " << emergenceAngle << std::endl;
}

  

  //calculate where the shower max is, which is from the ice in the direction of the shower pathLength away
  TVector3 xMaxPosition = destination+cosmicRayDirection*pathLength;


  if (debug) {
    std::cout << "anitaPosition: (" << anitaPosition.X() << "," << anitaPosition.Y() << "," << anitaPosition.Z() << ")" << std::endl;
    std::cout << "destination: (" << destination.X() << "," << destination.Y() << "," << destination.Z() << ")" << std::endl;
    std::cout << "destinationToSource: (" << destinationToSource.X() << "," << destinationToSource.Y() << "," << destinationToSource.Z() << ")" << std::endl;
    std::cout << "cosmicRayDirection: (" << cosmicRayDirection.X() << "," << cosmicRayDirection.Y() << "," << cosmicRayDirection.Z() << ")" << std::endl;
    std::cout << "xMaxPosition: (" << xMaxPosition.X() << "," << xMaxPosition.Y() << "," << xMaxPosition.Z() << ")" << std::endl;
  }


  // Calculate the geo-magnetic field field at x-max
  FieldPoint fp(usefulPat.realTime, xMaxPosition);
  if (debug) {
    std::cout << "FieldPoint: pos=(" << fp.posX() << "," << fp.posY() << "," << fp.posZ() << ")" << std::endl;
    std::cout << "           comp=(" << fp.componentX() << "," << fp.componentY() << "," << fp.componentZ() << std::endl;
  }


  // This is our electric field vector!
  // If I cared about getting the magnitude correct in addition to the orientation, there are some missing factors
  // It should be: B_vec x S_vec = (1/mu0)B^{2} E_vec
  // But there's no radio cherenkov/geomagnetic shower model or anything so for now just the cross product will do.
  TVector3 EVec = fp.field().Cross(cosmicRayDirection).Unit();
  if (debug) {
    std::cout << "EVec: (" << EVec.X() << "," << EVec.Y() << "," << EVec.Z() << ")" << std::endl;
  }
  
  // Here I find the VPol and HPol antenna axes.
  // And I'm going to  pretend that one of ANITA's antennas points exactly at phiWave
  // getting an off axis response will be more complicated

  // Since the antennas points down at -10 degrees, the VPol axis is 80 degrees above the horizontal plane
  TVector3 vPolAxis = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave, 80*TMath::DegToRad());
  // The VPol feed is up... (if) the HPol feed is to the right (looking down the boresight) then it points anticlockwise around the payload
  // phi increases anti-clockwise in payload coordinates, therefore
  TVector3 hPolAxis = getUnitVectorAlongThetaWavePhiWave(usefulPat, phiWave + TMath::PiOver2(), 0);
  
  // Dot the electric field with the antenna polarisation vectors...
  double vPolComponent = EVec.Dot(vPolAxis);
  double hPolComponent = EVec.Dot(hPolAxis);
  if (debug) {
    std::cout << "vPolComponent: " << vPolComponent << std::endl;
    std::cout << "hPolComponent: " << hPolComponent << std::endl;
  }
  
  // et voila
  double polarisationAngle = TMath::ATan(vPolComponent/hPolComponent);
  
  return polarisationAngle;


}




/** 
 * Plot the atmospheric density as a function of altitude used in the model
 * 
 * @return the canvas on which the plot is drawn
 */

TCanvas* GeoMagnetic::plotAtmosphere(){
  prepareGeoMagnetics();

  static int maxCanvas = 0;
  TCanvas* c1 = NULL;
  if(maxCanvas < 10){
    c1 = new TCanvas();
    grAtmosDensity.Draw("al");
    maxCanvas++;
  }
  return c1;
}



double GeoMagnetic::getAtmosphericDensity(double altitude){
  prepareGeoMagnetics();
  // TODO, convert altitude (above geoid) to height above MSL...
  return fExpAtmos->Eval(altitude);
}





/** 
 * Takes the final position and incoming direction and solves for the top of the atmosphere 
 * Currently uses 80km as the top of the atmosphere.
 * 
 * @param destination is the final position of the ray (either ANITA or the reflection point)
 * @param destinationToSource point away from the final position towards the arrival direction
 * 
 * @return the cartesian coordinate position vector at the top of the atmosphere
 */
TVector3 GeoMagnetic::getInitialPosition(const TVector3& destination, const TVector3& destinationToSource){
  prepareGeoMagnetics();

  if(TMath::Abs(destinationToSource.Mag() -1) > 1e-14){
    std::cerr  << "Warning in " << __PRETTY_FUNCTION__ << ", was expecting a unit vector, didn't get one. "
               << "Got " << 1.0-destinationToSource.Mag() << " away from mag==1... This calculation might be wonky... " << std::endl;
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





/** 
 * Takes a vector from the top of the atmosphere, and the direction of travel of the cosmic ray
 * and moves along the vector, integrating the atmospheric density until X_{max} is reached.
 * 
 * @param initialPosition is the cosmic ray's entry position to the atmosphere (cartesian TVector3)
 * @param cosmicRayDirection is it's direction of travel (cartesian TVector3)
 * 
 * @return a cartesian vector containing the position of X_{max}
 */
TVector3 GeoMagnetic::getXMaxPosition(const TVector3& initialPosition, const TVector3& cosmicRayDirection,
				      double xMax){
  prepareGeoMagnetics();

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
      std::cout << "Finding xMax position... " << xMax << " Initial position... " << currentLon << "\t" << currentLat << "\t" << currentAlt << "\t" << currentAtmosphereTraversed << std::endl;
    }
    
    currentAtmosphereTraversed += getAtmosphericDensity(currentAlt)*dx;

    if(debug && currentAtmosphereTraversed >= xMax){
      std::cout << "Found xMax position... " << xMax << " Initial position... " << currentLon << "\t" << currentLat << "\t" << currentAlt << "\t" << currentAtmosphereTraversed << std::endl;
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
    static int maxCanvas = 0;
    TCanvas* c1 = NULL;
    if(maxCanvas < 10){    
      c1 = new TCanvas();
      grAltPath->SetTitle("Cosmic Ray Altitude vs. distance traversed; Distance through atmosphere (m); Altitude (m)");
      grAltPath->SetBit(kCanDelete);
      grAltPath->Draw("al");
      maxCanvas++;
    }
  }
  
  return currentPosition;
}


