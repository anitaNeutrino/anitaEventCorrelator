#include "PayloadParameters.h"
#include "TMath.h"
#include "Adu5Pat.h"
#include "AntarcticaGeometry.h"
#include "FFTtools.h"

ClassImp(PayloadParameters); 

#ifdef USE_GEOGRAPHIC_LIB
#include <GeographicLib/GeodesicLine.hpp>
#endif

PayloadParameters::PayloadParameters(const Adu5Pat * gps, const AntarcticCoord & source_pos, const Refraction::Model * refract) 
  : payload(AntarcticCoord::WGS84, gps->latitude, gps->longitude, gps->altitude), 
    source(source_pos)
{

  payload.to(AntarcticCoord::CARTESIAN); 
  //vector from source to payload

  TVector3 p = payload.v(); 
  TVector3 s = source.v(); 

#ifndef USE_GEOGRAPHIC_LIB
  TVector3 sprime = s; 
  // this is stolen from AnitaGeomTool. Seems like it should be equivalent to a some form of rotateUz 
  // wihch would probably be a more efficient way to do this without trig functions 
  sprime.RotateZ(-1 * p.Phi()); 
  sprime.RotateY(-1 * p.Theta()); 
  sprime[2]=p.Mag()-fabs(sprime.z()); 

  sprime.RotateZ(gps->heading*TMath::DegToRad()); 

  //TODO: check if these axes need modification. Right now pitch and roll are 0 :) 
  sprime.RotateY(-gps->pitch *TMath::DegToRad()); 
  sprime.RotateX(-gps->roll *TMath::DegToRad()); 

  source_phi = FFTtools::wrap(sprime.Phi() * TMath::RadToDeg(),360); 
  source_theta = 90 - sprime.Theta() * TMath::RadToDeg(); 

  //Now do the opposite 
  TVector3 pprime = p;
  pprime.RotateZ(-1*s.Phi()); 
  pprime.RotateY(-1*s.Theta()); 
  pprime[2] = s.Mag() - fabs(pprime.z()); 
  payload_el = -FFTtools::wrap( 90 - pprime.Theta() * TMath::RadToDeg(), 180, 0); 
  payload_az = pprime.Phi() * TMath::RadToDeg(); 

#else

  //vector from source to paylaod
  TVector3 v = (s-p); 
  //angle between v and p gives zenith angle. elevation is - 90. 
  source_theta = p.Angle(v) * TMath::RadToDeg() - 90; 
  payload_el =  s.Angle(v) * TMath::RadToDeg() - 90; 
  //To get phi, we have to solve the inverse geodesic problem 
  AntarcticCoord swgs84 = source.as(AntarcticCoord::WGS84); 
  
  GeographicLib::Geodesic::WGS84().Inverse(gps->latitude, gps->longitude, swgs84.x, swgs84.y, source_phi, payload_az); 
  //rotate by 180 to get direction towards payload, and wrap 
  payload_az = FFTtools::wrap(payload_az-180, 360); 


  source_phi = FFTtools::wrap(gps->heading - source_phi,360); 
#endif

  if (refract) 
  {
    //We actually want the apparent anagle
    double payload_el_correction = 0;
    source_theta -= refract->getElevationCorrection(gps, &source, &payload_el_correction); 
    payload_el-= payload_el_correction; 
  }



  distance = (source.v() - payload.v()).Mag(); 
}



static const CartesianSurfaceMap & cartmap ( RampdemReader::dataSet d ) 
{
  if (d == RampdemReader::surface) 
  {
    static CartesianSurfaceMap sm(1000, d); 
    return sm; 
  }

  static CartesianSurfaceMap sm; 
  return sm; 
}



bool PayloadParameters::checkForCollision(double dx, AntarcticCoord * w, AntarcticCoord * w_exit,  RampdemReader::dataSet d, double grace, bool reverse) const
{

  AntarcticCoord x = (reverse ? payload: source).as(AntarcticCoord::CARTESIAN); 
  TVector3 v = (reverse ? -1 : 1) * (payload.v() - source.v()).Unit() * dx;  


  while(true)
  {
    x.to(AntarcticCoord::CARTESIAN); 
    x.x+= v.x();
    x.y+= v.y();
    x.z+= v.z(); 
 
    double height_above_surface = cartmap(d).metersAboveIce(x.x, x.y, x.z);

    if ( (!reverse && height_above_surface> 5000) || (reverse && (x.v() - source.v()).Mag2() < dx*dx))
      break; 

    if (height_above_surface  + grace < 0 )
    {
//      printf("BOOM! alt(%g,%g,%g)= %g\n", x.x, x.y, x.z, cartmap().surface(x.x, x.y)); 
      if (w) 
      {
        *w  = x; 
      }

      //figure out where no longer a collision
      if ( w_exit) 
      {
        while (true) 
        {
          x.x+= v.x();
          x.y+= v.y();
          x.z+= v.z(); 
          height_above_surface = cartmap(d).metersAboveIce(x.x, x.y, x.z);

          if (height_above_surface + grace >=0 ) 
          {
            *w_exit = x; 
            break; 
          }
          else if ( reverse && (x.v() - source.v()).Mag2() < dx*dx)
          {
            *w_exit = source; 
            break; 
          }

        }
      }

      return true; 
    }
  }


  return false; 


}


PayloadParameters::PayloadParameters(const PayloadParameters & other) 
{
  source_phi = other.source_phi; 
  source_theta = other.source_theta; 
  payload_el = other.payload_el; 
  payload_az = other.payload_az; 
  distance = other.distance; 
  payload = other.payload; 
  source = other.source; 
}



int PayloadParameters::findSourceOnContinent(double theta, double phi, const Adu5Pat * gps, PayloadParameters * p, 
                                             const Refraction::Model * m, 
                                             double collision_check_dx,
                                             double min_dx, double tol, double min_el, 
                                             RampdemReader::dataSet d) 
{
  //no chance. 
  if (theta < 0 ) 
  {
    return 0; 
  }

  AntarcticCoord payload(AntarcticCoord::WGS84, gps->latitude, gps->longitude, gps->altitude); 
  AntarcticCoord x = payload.as(AntarcticCoord::CARTESIAN); 


#ifdef USE_GEOGRAPHIC_LIB
  GeographicLib::GeodesicLine gl(GeographicLib::Geodesic::WGS84(),  gps->latitude,  gps->longitude, gps->heading - phi); 
  
  size_t i = 1; 
  double step = 1000*300;  //start with 300 km step 
  PayloadParameters ok;
  int loopCount = 0;
  while (loopCount<100) //simply hard coded the upper limit of the loop, otherwise it may stuck in the loop forever.
  {
    loopCount++;
    double lat, lon; 
    gl.Position(i * step, lat, lon); 
    AntarcticCoord c(AntarcticCoord::WGS84, lat, lon, RampdemReader::SurfaceAboveGeoid(lon,lat,d)); 
    *p = PayloadParameters(gps, c, m); 

    //printf("%f::%f %f::%f::  %g   %f\n", phi, p->source_phi, theta, p->source_theta, p->payload_el, i * step); 


    //we found something that works
    if (fabs(p->source_phi -phi) < tol && fabs(p->source_theta - theta) < tol && p->payload_el >= min_el) 
    {
      //check for a collision, 
      //then go to exit point? 
      if (collision_check_dx && p->checkForCollision( TMath::Min(collision_check_dx,step), 0, &c,   d))
      {
        step= TMath::Min(collision_check_dx, step) ; 
        double dist = AntarcticCoord::surfaceDistance(payload, c); 
        i = dist/step; 
        step = dist/i; 
        continue; 
      }

      return 1; 
    }

    //overshot the source, or are over the horizon with too large a step size , we should lower the step size and go back 
    if (p->source_theta < theta || (p->payload_el < min_el && step > min_dx)) 
    {
      step = step / 2; 
      i = i*2-1; 
      continue; 
    }
    else if (p->payload_el < min_el)  
    {
      *p = ok; 
      return 0; 
    }
    else if ( step < min_dx) 
    {
      return -1; 
    }

    //store last ok thing) 
    if (p->payload_el >= min_el) 
    {
      ok = *p; 
    }

    i++; 
  }



 #else

  /* This doesn't work yet 
  //use great ellipse 

  //Find vector normal to our plane, use 0,0,0 as our point
  TVector3 d = (x.x,x.y,0); //north pointing vector
  d.RotateZ(gps->heading); 

  TVector3 n = x.v().Cross(d); 

  // The intersection of the plane and geoid will be an ellipse. We must follow it. 
  // This sucks. 

  while(true) 
  {
    //somehow propagate along ellipse 

    //form the coordinates on the ice surface 
    AntarcticCoord s(AntarcticCoord::CARTESIAN, x2.X(), x2.Y(), cartmap().z(x2.X(), x2.Y()));
    TString str;
    s.asString(&str); 
    printf("%s\n",str.Data());

    PayloadParameters pp(x,gps->heading, s); 

    printf("%f::%f %f::%f\n", phi, pp.source_phi, theta, pp.source_theta); 
    if (isnan(s.z)) return 0; 

    if (fabs(pp.source_phi -phi) < tol && fabs(pp.source_theta - theta) < tol) 
    {
      return new PayloadParameters(pp); 
    }
  }
  */
#endif

  return 0; 
}
