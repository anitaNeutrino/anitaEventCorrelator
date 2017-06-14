#include "GeoMagnetic.h"
#include "TH2D.h"
#include "AnitaGeomTool.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "AntarcticaBackground.h"


void testIGRF(){

  GeoMagnetic::setDebug(true);
  
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
      double Z_lla = GeoMagnetic::Z_atLonLatAlt(0, lon, lat, 40e3);
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
