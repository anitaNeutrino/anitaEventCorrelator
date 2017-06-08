#include "GeoMagnetic.h"
#include "TH2D.h"
#include "AnitaGeomTool.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "AntarcticaBackground.h"

void lonLatAltToSpherical(double lon, double lat, double alt, double& r, double& theta, double& phi){
  double cartesian[3];
  AnitaGeomTool::Instance()->getCartesianCoords(lat, lon, alt, cartesian);
  double x = cartesian[0];
  double y = cartesian[1];
  double z = cartesian[2];

  // AnitaGeomTool has something really sinster going on?
  z = lat >= 0 ? TMath::Abs(z) : -TMath::Abs(z);

  r = TMath::Sqrt(x*x + y*y + z*z);
  theta = r > 0 ? TMath::ACos(z/r) : 0;
  phi = -TMath::ATan2(y, x) + 0.5*TMath::Pi();
  phi = phi >= TMath::Pi() ?  phi - TMath::TwoPi() : phi;

  // std::cout << lon <<  "\t" << lat << "\t" << alt  << "\t" <<  r  << "\t" << theta << "\t" << phi << std::endl;
  // std::cout << lon <<  "\t" << lat << "\t" << alt  << "\t" <<  r  << "\t" << theta << "\t" << phi << "\t" <<  z << std::endl;  
}

void testIGRF(){

  const int nx = 360;
  const int ny = 180;
  auto hS = new TH2D("hS", "Z component of geo-magnetic field (spherical coordinates)",
                     nx, -180, 180,
                     ny, 0,  180);
                    // ny, -180,  0);
  auto hL = new TH2D("hL", "Z component of geo-magnetic field (lon/lat/alt)",
                     nx, -180, 180,
                     ny, -90,  90);
                    // ny, -180,  0);
  
  const double r = 6371.2e3;
    
  for(int by=1; by <= ny;  by++){
    double lat = hL->GetYaxis()->GetBinLowEdge(by);
    double theta = hS->GetYaxis()->GetBinLowEdge(by)*TMath::DegToRad();
    for(int bx=1; bx <= nx;  bx++){
      double lon = hL->GetXaxis()->GetBinLowEdge(bx);
      double phi = hS->GetXaxis()->GetBinLowEdge(bx)*TMath::DegToRad();
      double Z_spherical = GeoMagnetic::Z_atSpherical(0, r, theta, phi);
      hS->SetBinContent(bx, by, Z_spherical);
      double Z_lla = GeoMagnetic::Z_atLonLatAlt(0, lon, lat, 0);
      std::cout << Z_lla << std::endl;
      hL->SetBinContent(bx, by, Z_lla);
      
    }
  }
  auto c1 = new TCanvas();
  c1->Divide(2);
  c1->cd(1);
  hS->Draw("colz");
  c1->cd(2);
  hL->Draw("colz");

  auto c2 = new TCanvas();  
  TGraph* grTheta = new TGraph();
  for(int lat=-90; lat <=90; lat++){
    double lon = 0, alt = 0, r, theta, phi;
    lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
    grTheta->SetPoint(grTheta->GetN(), lat, theta*TMath::RadToDeg());
  }
  grTheta->SetTitle("Latitude tranform; Latitude (Deg); #theta (Deg)");
  grTheta->Draw();

  auto c3 = new TCanvas();  
  TGraph* grPhi = new TGraph();
  for(int lon=-180; lon <=180; lon++){
    double lat = 0, alt = 0, r, theta, phi;
    lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
    grPhi->SetPoint(grPhi->GetN(), lon, phi*TMath::RadToDeg());
  }
  grPhi->SetTitle("Longitude tranform; Longitude (Deg); #phi (Deg)");
  grPhi->Draw();

  
  // auto c4 = new TCanvas();
  // auto x = new AntarcticaBackground();
  // x->Draw();
}

















