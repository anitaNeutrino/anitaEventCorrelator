

void testAstro()
{


  AnitaDataset d(342); 
  UsefulAdu5Pat pat(d.gps()); 

  double dec, ra; 



  for (double phi = 0; phi < 360; phi++)
  {
    for (double theta = -30; theta < 30; theta++)
    {
      pat.astronomicalCoordinates(phi,theta,&ra,&dec); 
      double p,t; 
      pat.fromRADec(ra,dec,&p,&t); 
      printf("%g %g -> %g %g -> %g %g\n",phi,theta,ra,dec,p,t); 

    }

  }


}
