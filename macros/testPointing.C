void testPointing(int run = 342, int event = 58023120) 
{

  AnitaDataset d(run); 
  d.getEvent(event); 

  UsefulAdu5Pat gps(d.gps()); 

  double wais_phi, wais_theta; 
  gps.getThetaAndPhiWaveWaisDivide(wais_theta, wais_phi); 

  printf("Real WAIS lat/lon/alt: %g %g %g\n", AnitaLocations::getWaisLatitude(), AnitaLocations::getWaisLongitude(), AnitaLocations::getWaisAltitude()); 
  printf("WAIS at phi=%g theta=%g\n", wais_phi* TMath::RadToDeg(), wais_theta* TMath::RadToDeg()); 


  double lat,lon,alt; 
  int trace = gps.traceBackToContinent(wais_phi, wais_theta, &lon, &lat,&alt,0); 

  printf("traceBackToContent projects to lat/lon/alt: %g %g %g\n",lat,lon,alt); 



  double payload_phi, payload_theta; 
  gps.getThetaAndPhiWave(lon,lat,alt, payload_theta, payload_phi); 
  printf("In payload coordinates that is: phi=%g theta=%g\n", payload_phi * TMath::RadToDeg(), payload_theta * TMath::RadToDeg()); 


  AntarcticCoord c(AntarcticCoord::WGS84,lat,lon,alt); 
  PayloadParameters pp(d.gps(), c); 
  printf("Payload parameters thinks it's at phi=%g, theta=%g\n", pp.source_phi, pp.source_theta); 

  int get = gps.getSourceLonAndLatAtAlt(wais_phi, wais_theta, lon, lat, alt);
  printf("getSourceLonAndLatAtAlt projects to lat/lon/alt: %g %g %g\n",lat,lon,alt); 


  gps.getThetaAndPhiWave(lon,lat,alt, payload_theta, payload_phi); 
  printf("In payload coordinates that is: phi=%g theta=%g\n", payload_phi * TMath::RadToDeg(), payload_theta * TMath::RadToDeg()); 

  AntarcticCoord c2(AntarcticCoord::WGS84,lat,lon,alt); 
  PayloadParameters pp2(d.gps(), c2); 
  printf("Payload parameters thinks it's at phi=%g, theta=%g\n", pp2.source_phi, pp2.source_theta); 



}
