

void testConversion()
{

  for (double lat = -90; lat <= -60; lat+=5) 
  {
    for (double lon = -180; lon <= 180; lon+=5) 
    {
      printf("Lat: %g, Lon: %g\n", lat,lon); 

      AntarcticCoord c(AntarcticCoord::WGS84, lat, lon, 1000); 

      AntarcticCoord direct = c.as(AntarcticCoord::STEREOGRAPHIC); 
      printf("  Direct Stereographic: (%g %g %g)\n", direct.x, direct.y, direct.z); 

      AntarcticCoord indirect = c.as(AntarcticCoord::CARTESIAN).as(AntarcticCoord::STEREOGRAPHIC); 
      printf("  Indirect Stereographic: (%g %g %g)\n", indirect.x, indirect.y, indirect.z); 

      printf("  Delta:  (%g, %g, %g)\n", direct.x - indirect.x, direct.y - indirect.y, direct.z - indirect.z); 
    }
  }



}
