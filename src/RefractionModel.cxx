#include "RefractionModel.h" 
#include "AntarcticAtmosphere.h" 
#include "AntarcticaGeometry.h" 
#include "TEllipse.h" 
#include "Adu5Pat.h" 
#include "TH2.h" 
#include "TGraph.h" 


const int nH = 5; 
const int nh = 4; 


static const double A[nh][nH] = 
{ 
  { 3.5553, 4.7295, 4.4731, 5.1151, 4.8049 },
  { 3.2987, 3.1228, 4.0890, 3.7471, 4.4854 },
  { 2.9440, 2.8858, 2.6563, 3.4098, 3.2843 },
  { 2.4514, 2.6189, 2.4666, 2.8971, 3.1265 }
};


static const double B[nh][nH] = 
{ 
  { -0.0095, 0.1173, 0.0907, 0.1432, 0.1149 },
  { 0.0029, -0.0141, 0.0898, 0.0533, 0.1190 },
  { 0.0001, -0.0052, -0.0245, 0.0517, 0.0406 },
  { -0.0156, -0.0026, -0.0170, 0.0214, 0.0584 }
};


static const double C[nh][nH] = 
{
  { 4.3780, 4.5653, 4.5795, 4.7217, 4.7491 }, 
  { 4.2998, 4.3507, 4.5034, 4.5180, 4.6733 },
  { 4.2362, 4.2858, 4.3466, 4.4356, 4.4903 },
  { 4.3714, 4.3139, 4.3234, 4.3701, 4.4639 }
};

static const double D[nh][nH] = 
{
  { 5.0281, 7.0018, 6.3010, 7.9669, 7.1021 }, 
  { 5.4476, 5.3225, 6.8895, 5.8489, 7.8475 },
  { 5.5789, 5.8992, 5.0228, 6.4541, 6.1882 },
  { 6.4298, 7.4389, 7.0311, 7.2302, 7.9347 }
};

static const double E[nh][nH] = 
{
  { -5.5423, -5.8028, -5.8366, -6.0156, -6.0501 }, 
  { -5.4901, -5.5646, -5.7826, -5.7774, -5.9982 }, 
  { -5.4010, -5.5246, -5.5159, -5.7338, -5.7958 }, 
  { -5.3860, -5.5437, -5.6085, -5.7002, -5.8241 } 
};

static TH2 * hCoeffs[5] = {0}; 

__attribute__ ((constructor)) 
static void init_hists()
{
  for (int coeff = 0; coeff < 5; coeff++)
  {

     hCoeffs[coeff] = new TH2D(TString::Format("pg_coeffs_%c", 'A' + coeff),TString::Format("%c Coeffs", 'A' + coeff) , 
                            nh, -0.5e3, 0.35e3, nH, 35.5e3 , 40.5e3); 

     hCoeffs[coeff]->SetDirectory(0); 

     for (int i = 1; i <= nh; i++) 
     {
       for (int j = 1; j<= nH; j++)
       {
             hCoeffs[coeff]->SetBinContent(i,j,  (coeff == 0 ? A :
                                                 coeff == 1 ? B :
                                                 coeff == 2 ? C :
                                                 coeff == 3 ? D :
                                                 E)[i-1][j-1]); 
       }
    }
  }
}



double Refraction::PositionIndependentModel::getElevationCorrection(const Adu5Pat * pat, const AntarcticCoord * source, double * correction_at_source)  const
{

  PayloadParameters pp(pat,*source); 
  double H = AntarcticAtmosphere::WGS84toMSL(pat); 
  double h = AntarcticAtmosphere::WGS84toMSL(pat); 
  return getElevationCorrection(pp.source_theta, h,H, correction_at_source); 

}

double Refraction::PGFit::getElevationCorrection(double theta, double h, double H, double * correction_at_source)  const
{

  // pin to limits
   if (h > 4e3) h = 4e3; 
   if (h < 0) h = 0; 
   if (H < 36e3) H = 36e3; 
   if (H > 40e3) H = 40e3; 
  
   double a = hCoeffs[0]->Interpolate(h,H); 
   double b = hCoeffs[1]->Interpolate(h,H); 
   double c = hCoeffs[2]->Interpolate(h,H); 
   double d = hCoeffs[3]->Interpolate(h,H); 
   double e = hCoeffs[4]->Interpolate(h,H); 
   double el = -theta; 
   return -(a / (el * el) + b / (el) + c * exp(d * (el-e))); 
   if (correction_at_source) *correction_at_source = 0; //we don't provide one 
} 


const double REARTH = 6371e3; 
const double speed_of_light=299792458; 
int Refraction::RaytracerSpherical::raytrace(const Setup * setup, Result * result) 
{

  last_path_x.clear(); 
  last_path_y.clear(); 

  double r0 = setup->start_alt + REARTH; 
  double r1 = setup->end_alt + REARTH; 
  double theta = setup->thrown_payload_angle * TMath::DegToRad(); 
  double r = r0; 
  double phi=0; 
  double x0 = 0; 
  double y0 = r0; 

  double x= x0; 
  double y= y0; 
  result->ray_distance = 0; 
  result->ray_time = 0; 

  while (r <= r1) 
  {

    last_path_x.push_back(x); 
    last_path_y.push_back(y); 

    // we want to step in r  dphi 

    double dphi = step_size  / r; 
    double dr= step_size * tan(theta); 

    double N = m->get(r-REARTH, AntarcticAtmosphere::REFRACTIVITY); 
    double N1 = m->get(r-REARTH+dr, AntarcticAtmosphere::REFRACTIVITY); 
    double DN = N1-N ;
    double dn = DN * 1e-6; 
    double n = 1 + N * 1e-6; 
 //   printf("%g\n",tan(theta) * r); 
    double dtheta = dn/n/tan(theta) + dphi; 
//    printf("%g %g %g\n", r,theta,dtheta); 

    double dL = sqrt((r*dphi)*(r*dphi) + dr*dr); 
    result->ray_distance += dL; 
    result->ray_time += dL * n / speed_of_light; 
//    printf("%g %g %g %g\n",dphi,dtheta, phi,theta); 
    theta += dtheta; 
    phi +=  dphi; 
    r += dr; 
    x = r* sin(phi); 
    y = r* cos(phi); 
  }

  last_path_x.push_back(x); 
  last_path_y.push_back(y); 

  double dx = x- x0; 
  double dy = y- y0; 
  double dM  = sqrt(dx*dx + dy*dy); 

  result->actual_distance = dM; 
  result->apparent_source_angle = theta * TMath::RadToDeg(); 
  result->actual_source_angle = 90 - acos( (dx*x + dy*y) / (r * dM)) *TMath::RadToDeg(); 
  result->actual_payload_angle = 90 - acos(dy/dM) * TMath::RadToDeg(); 
  return 0; 
}



void Refraction::RaytracerSpherical::draw() 
{
  TGraph * g = new TGraph (last_path_x.size(), &last_path_x[0], &last_path_y[0]); 
  g->Draw("alp"); 

  //TEllipse * e = new TEllipse(0,0,REARTH, REARTH); 
//  e->SetLineColor(4); 
//  e->Draw("lsame"); 

}



double Refraction::SphRay::getElevationCorrection(double el, double hSource, double hPayload, double * correction_at_source ) const
{


  UInt_t hash = 0;

  if (use_cache) 
  {
    //use a very dumb hash for now, assuming a certain range for the start and end heights and 10 m  / 0.01 degree height / angle precision. 

    //box hSource into (0,5000,nearest multiple of 100) 
    hSource = TMath::Max(0., hSource); 
    hSource = TMath::Min(4999., hSource); 
    hSource = ((int)hSource/100) * 100; 
    hash += hSource / 100; 

     //box hPayload into (36000,41000,nearest multiple of 100) 
    hPayload = TMath::Max(36e3, hPayload); 
    hPayload = TMath::Min(40.999e3, hPayload); 
    hPayload = ((int)hPayload/100) * 100; 
    hash +=  50 * ((hPayload-36e3) / 100); 

    
    // box el into (0,90, nearest hundredth of a degree) 
    el = TMath::Min(90.,el); 
    el = TMath::Max(0.,el); 
    el = ((int) (el * 100)) / 100.; 
//    printf("%g %g %g\n", hSource, hPayload, el); 
    hash += 50 * 50 * el * 100; 

    if (cache.count(hash)) 
    {
      if (correction_at_source) 
        *correction_at_source = cache[hash].second; 
      return cache[hash].first; 
    }
  }

  //create the ray tracer 
  RaytracerSpherical ray(atm); 
  ray.step_size = step;  

  //compute the horizon angle at the ground 
  
  double throw_angle = TMath::RadToDeg() * acos (  cos(el * TMath::DegToRad())  * (REARTH + hPayload) / (REARTH + hSource)  * (1 + 1e-6 * atm->get(hPayload, AntarcticAtmosphere::REFRACTIVITY)) /  (1 + 1e-6 * atm->get(hSource, AntarcticAtmosphere::REFRACTIVITY)));

  RaytracerSpherical::Setup s; 
  RaytracerSpherical::Result r; 
  s.save_path = false; 
  s.start_alt = hSource; 
  s.end_alt = hPayload; 
  s.thrown_payload_angle = throw_angle; 

  ray.raytrace(&s,&r); 


  double answer = r.actual_source_angle - r.apparent_source_angle; 
  double correction = r.actual_payload_angle - s.thrown_payload_angle; 
  

  if (use_cache) 
  {
    cache[hash] = std::pair<double,double> (answer, correction) ; 
  }

  if (correction_at_source) *correction_at_source = correction; 
  return answer ; 
}







