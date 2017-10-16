#include "RefractionModel.h" 

void testInversion(double hmin = 0, double hmax = 200, double payload=35e3, double A =0.1) 
{
  
  const AntarcticAtmosphere::AtmosphericModel & m = AntarcticAtmosphere::ITURefractivity(); 
  const AntarcticAtmosphere::ArtificialInversion & im = AntarcticAtmosphere::ArtificialInversion(m, hmax, A); 
  TCanvas * c = new TCanvas; 


  c->Divide(2,1); 
  c->cd(1); 
  im.makeGraph(hmin,payload, 1000, AntarcticAtmosphere::REFRACTIVITY)->Draw("alp"); 

  Refraction::RaytracerSpherical ray(&m); 

  ray.step_size=10; 
  Refraction::RaytracerSpherical::Setup * s =  new Refraction::RaytracerSpherical::Setup;; 
  Refraction::RaytracerSpherical::Result * r  = new Refraction::RaytracerSpherical::Result; 

  TFile f("inversion.root","RECREATE"); 
  TTree * t = new TTree("raytrace","Raytrace"); 
  t->Branch("setup",&s); 
  t->Branch("result",&r); 
  c->cd(2); 

  double h = hmin; 
  double H = payload; 
  double step = 0.0005; 
  for (double el = 0.0005; el <= 89; el+= step)
  {
        if (el >= 0.1) step= 0.005; 
        if (el >= 1)  step = 0.05;

        s->start_alt = h; 
        s->end_alt = H; 
        s->thrown_payload_angle = el;
        ray.raytrace(s,r); 
        t->Fill(); 


        if (el == 0.0005) ray.draw(); 


  }

  t->Write(); 


}


