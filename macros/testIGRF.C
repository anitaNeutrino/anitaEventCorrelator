#include "GeoMagnetic.h"
#include "TH2D.h"
#include "AnitaGeomTool.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "AntarcticaBackground.h"
#include "TDatime.h"



/** 
 * Convert unixTime to year
 * 
 * @param unixTime is the number of seconds since 1970
 * 
 * @return year as a decimal quantity
 */
std::pair<int, double> unixTimeToFractionalYear2(UInt_t unixTime){

  std::map<UInt_t, Int_t> unixTimeToYear; // map of unixTime to year (1st Jan zero seconds past midnight GMT), every fifth year
  unixTimeToYear[0] = 1970;
  unixTimeToYear[157766400] = 1975;
  unixTimeToYear[315532800] = 1980;
  unixTimeToYear[631152000] = 1990;
  unixTimeToYear[788918400] = 1995;
  unixTimeToYear[946684800] = 2000;
  unixTimeToYear[1104537600] = 2005;
  unixTimeToYear[1262304000] = 2010;
  unixTimeToYear[1420070400] = 2015;
  unixTimeToYear[1577836800] = 2020;
  unixTimeToYear[1735689600] = 2025;
  unixTimeToYear[1893456000] = 2030;
  unixTimeToYear[2051222400] = 2035;
  
  std::map<UInt_t, Int_t>::iterator next5YearIt = unixTimeToYear.upper_bound(unixTime);

  std::pair<int, double> fifthYearAndFrac;

  if(next5YearIt != unixTimeToYear.begin()){

    // upper bound returns first item greater than the search
    UInt_t nextUnixTime = next5YearIt->first;

    // we also want the last time less than or equal to, so decrement iterator
    std::map<UInt_t, Int_t>::iterator last5YearIt = next5YearIt;
    last5YearIt--;

    // upper bound returns first item greater than the search
    Int_t last5Year = last5YearIt->second;
    UInt_t lastUnixTime = last5YearIt->first;

    double fracThrough5Year = double(unixTime - lastUnixTime)/double(nextUnixTime - lastUnixTime);
    
    return std::pair<int, double>(last5Year, fracThrough5Year);
    
  }
  else{
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", the laziness of a programmer in the distant past has caused you a problem... could not figure out the year from unixTime!" << std::endl;
    return std::pair<int, double> (0, 0);
  }  
}



void testIGRF(){


  std::vector<UInt_t> unixTimes {1418869215, 1480520956, 1482518716, 1480687576, 1420170374, 1893456000-1};  

  for(auto& unixTime : unixTimes){

    std::pair<int, double> fracYear2 = unixTimeToFractionalYear2(unixTime);    
    std::cout << fracYear2.first << "\t" << fracYear2.second << std::endl;
  }
  
  
  // return;

  // GeoMagnetic::setDebug(true);
  
  GeoMagnetic::plotAtmosphere();

  Adu5Pat pat;
  pat.longitude = 139;
  pat.altitude =  40e3;
  pat.latitude = -79;

  pat.heading = 27;
  pat.pitch = 0;
  pat.roll = 0;

  UsefulAdu5Pat usefulPat(&pat);

  // double thetaWave = 1.11111 * TMath::DegToRad();
  double thetaWave = -2 * TMath::DegToRad();
  double phiWave = 13.333 * TMath::DegToRad();

  std::cout << "Trying to get phiWave = " << phiWave*TMath::RadToDeg() << ", thetaWave = " << thetaWave*TMath::RadToDeg() << std::endl;
  
  double polAngle = GeoMagnetic::getExpectedPolarisation(usefulPat, phiWave, thetaWave);
  std::cout << TMath::RadToDeg()*polAngle << std::endl;

  auto cc = GeoMagnetic::plotFieldAtAltitude(0, 40e3);
  
  const int nx = 360;
  const int ny = 180;
  auto hZ = new TH2D("hZ", "z-component of geo-magnetic field (spherical coordinates)",
                     nx, -180, 180,
                     ny, -90,  90);
  auto hY = new TH2D("hY", "x-component of geo-magnetic field (spherical coordinates)",
                     nx, -180, 180,
                     ny, -90,  90);
  auto hX = new TH2D("hX", "x-component of geo-magnetic field (spherical coordinates)",
                     nx, -180, 180,
                     ny, -90,  90);
  auto hL = new TH2D("hL", "Z component of geo-magnetic field (lon/lat/alt)",
                     nx, -180, 180,
                     ny, -90,  90);
  const double alt = 0;
  
  for(int by=1; by <= ny;  by++){
    double lat = hL->GetYaxis()->GetBinLowEdge(by);
    // double theta = hZ->GetYaxis()->GetBinLowEdge(by)*TMath::DegToRad();
    for(int bx=1; bx <= nx;  bx++){
      double lon = hL->GetXaxis()->GetBinLowEdge(bx);
      // double phi = hZ->GetXaxis()->GetBinLowEdge(bx)*TMath::DegToRad();      

      GeoMagnetic::FieldPoint f(0, lon, lat, alt);

      // hZ->SetBinContent(bx, by, f.componentZ());
      // hY->SetBinContent(bx, by, f.componentY());
      // hX->SetBinContent(bx, by, f.componentX());
      
      // hZ->SetBinContent(bx, by, f.fField.Mag());
      // hY->SetBinContent(bx, by, f.fField.Phi());
      // hX->SetBinContent(bx, by, f.fField.Theta());
      double Z_lla = GeoMagnetic::Z_atLonLatAlt(1893456000, lon, lat, 40e3);
      // std::cout << Z_lla << std::endl;
      hL->SetBinContent(bx, by, Z_lla);
    }
  }
  
  auto c1 = new TCanvas();
  hL->Draw("colz");

  // auto c1b = new TCanvas();
  // c1b->Divide(3);
  // c1b->cd(1);
  // hX->Draw("colz");
  // c1b->cd(2);
  // hY->Draw("colz");
  // c1b->cd(3);
  // hZ->Draw("colz");

  // return;

  // auto c2 = new TCanvas();  
  // TGraph* grTheta = new TGraph();
  // for(int lat=-90; lat <=90; lat++){
  //   double lon = 0, alt = 0, r, theta, phi;
  //   lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  //   grTheta->SetPoint(grTheta->GetN(), lat, theta*TMath::RadToDeg());
  // }
  // grTheta->SetTitle("Latitude tranform; Latitude (Deg); #theta (Deg)");
  // grTheta->Draw();

  // auto c3 = new TCanvas();  
  // TGraph* grPhi = new TGraph();
  // for(int lon=-180; lon <=180; lon++){
  //   double lat = 0, alt = 0, r, theta, phi;
  //   lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
  //   grPhi->SetPoint(grPhi->GetN(), lon, phi*TMath::RadToDeg());
  // }
  // grPhi->SetTitle("Longitude tranform; Longitude (Deg); #phi (Deg)");
  // grPhi->Draw();
  
            
}
