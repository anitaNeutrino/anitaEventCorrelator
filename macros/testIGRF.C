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

  // std::cout << lon << "\t" << lat << "\t" << alt << "\t" << x << "\t" << y << "\t" << z << std::endl;

  // AnitaGeomTool has something really sinster going on?
  z = lat >= 0 ? TMath::Abs(z) : -TMath::Abs(z);
  // z = -z;

  r = TMath::Sqrt(x*x + y*y + z*z);
  theta = r > 0 ? TMath::ACos(z/r) : 0;
  phi = -TMath::ATan2(y, x) + 0.5*TMath::Pi();
  phi = phi >= TMath::Pi() ? phi - TMath::TwoPi() : phi;

  // std::cout << lon <<  "\t" << lat << "\t" << alt  << "\t" <<  r  << "\t" << theta << "\t" << phi << std::endl;
  // std::cout << lon <<  "\t" << lat << "\t" << alt  << "\t" <<  r  << "\t" << theta << "\t" << phi << "\t" <<  z << std::endl;  
}



void sphericalToLatLonAlt(double& lon, double& lat, double& alt, double r, double theta, double phi){

  double x = r*TMath::Sin(phi)*TMath::Sin(theta);
  double y = r*TMath::Cos(phi)*TMath::Sin(theta);
  double z = r*TMath::Cos(theta);
  double cartesian[3] = {x, y, z};

  auto g = AnitaGeomTool::Instance();
  g->getLatLonAltFromCartesian(cartesian, lat, lon, alt);

  // fml
  lat = theta*TMath::RadToDeg() <= 90 ? -lat : lat;
}




void testIGRF(){
  {
    double lat=-30, lon=45, alt=40e3;
    double r, theta, phi;
    lonLatAltToSpherical(lon, lat, alt, r, theta, phi);
    std::cout << "input lon/lat/alt = " << lon << "/" << lat << "/" << alt << std::endl;
    std::cout << "output r/theta/phi = " << r << "/" << theta*TMath::RadToDeg() << "/" << phi*TMath::RadToDeg() << std::endl;
    sphericalToLatLonAlt(lon, lat, alt, r, theta, phi);
    std::cout << "output lon/lat/alt = " << lon << "/" << lat << "/" << alt << std::endl;
  }

  const int nx = 360;
  const int ny = 180;
  auto hZ = new TH2D("hZ", "Z component of geo-magnetic field (spherical coordinates)",
                     nx, -180, 180,
                     ny, -90,  90);
  auto hY = new TH2D("hY", "Y component of geo-magnetic field (spherical coordinates)",
                     nx, -180, 180,
                     ny, -90,  90);
  auto hX = new TH2D("hX", "X component of geo-magnetic field (spherical coordinates)",
                     nx, -180, 180,
                     ny, -90,  90);
  auto hL = new TH2D("hL", "Z component of geo-magnetic field (lon/lat/alt)",
                     nx, -180, 180,
                     ny, -90,  90);
  
  const double r = 6371.2e3;
  const double alt = 0;
  
  for(int by=1; by <= ny;  by++){
    double lat = hL->GetYaxis()->GetBinLowEdge(by);
    // double theta = hZ->GetYaxis()->GetBinLowEdge(by)*TMath::DegToRad();
    for(int bx=1; bx <= nx;  bx++){
      double lon = hL->GetXaxis()->GetBinLowEdge(bx);
      // double phi = hZ->GetXaxis()->GetBinLowEdge(bx)*TMath::DegToRad();      

      GeoMagnetic::Field f = GeoMagnetic::getFieldAtLonLatAlt(0, lon, lat, alt);
      
      hZ->SetBinContent(bx, by, f.Z);
      hY->SetBinContent(bx, by, f.Y);
      hX->SetBinContent(bx, by, f.X);
      double Z_lla = GeoMagnetic::Z_atLonLatAlt(0, lon, lat, 0);
      // std::cout << Z_lla << std::endl;
      hL->SetBinContent(bx, by, Z_lla);
    }
  }
  auto c1 = new TCanvas();
  hL->Draw("colz");

  auto c1b = new TCanvas();
  c1b->Divide(3);
  c1b->cd(1);
  hX->Draw("colz");
  c1b->cd(2);
  hY->Draw("colz");
  c1b->cd(3);
  hZ->Draw("colz");

  return;

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

  
  auto c4 = new TCanvas();
  auto x = new AntarcticaBackground();

  auto z = new TH2D("hVert", "",
                    x->GetNbinsX(), x->GetXaxis()->GetBinLowEdge(1), x->GetXaxis()->GetBinUpEdge(x->GetNbinsX()),
                    x->GetNbinsY(), x->GetYaxis()->GetBinLowEdge(1), x->GetYaxis()->GetBinUpEdge(x->GetNbinsY()));
  x->Draw();
}

















