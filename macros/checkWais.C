#include "AnitaDataset.h" 
#include "AntarcticaGeometry.h" 
#include "UsefulAdu5Pat.h" 
#include "AnitaConventions.h" 

void checkWais(int run = 342, int entry = 0) 
{

  AnitaDataset d(run); 
  d.getEntry(entry); 

  UsefulAdu5Pat gps(d.gps()); 

  double phi,theta; 
  gps.getThetaAndPhiWaveWaisDivide(theta,phi) ; 
  theta *= TMath::RadToDeg(); 
  phi *= TMath::RadToDeg(); 


  printf("UsefulAdu5Pat theta and phi: %f %f\n", theta, phi); 


  AntarcticCoord wais(AntarcticCoord::WGS84, AnitaLocations::LATITUDE_WAIS_A3, 
                                             AnitaLocations::LONGITUDE_WAIS_A3, 
                                             AnitaLocations::ALTITUDE_WAIS_A3);  

  PayloadParameters p(d.gps(), wais); 

  printf("Payload parameters source theta and phi: %f %f\n", p.source_theta, p.source_phi); 
  printf("Payload parameters distance, payload az,el: %f %f %f\n", p.distance, p.payload_az, p.payload_el); 


}
