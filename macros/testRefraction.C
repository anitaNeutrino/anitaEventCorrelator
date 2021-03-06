#include "RefractionModel.h" 

void testRefraction() 
{
  
  const AntarcticAtmosphere::AtmosphericModel & m = AntarcticAtmosphere::ITURefractivity(); 

  Refraction::RaytracerSpherical ray(&m); 

  ray.step_size=10; 
  Refraction::RaytracerSpherical::Setup * s =  new Refraction::RaytracerSpherical::Setup;; 
  Refraction::RaytracerSpherical::Result * r  = new Refraction::RaytracerSpherical::Result; 

  TFile f("raytrace.root","RECREATE"); 
  TTree * t = new TTree("raytrace","Raytrace"); 
  t->Branch("setup",&s); 
  t->Branch("result",&r); 

  for (double h = 0; h <= 4000; h+= 500)
  {
    for (double H = 35e3; H <= 40e3; H+= 500)
    {
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


        if (h == 0 && H == 40e3 && el == 0.0005) ray.draw(); 


      }
    }
  }

  t->Write(); 


}


