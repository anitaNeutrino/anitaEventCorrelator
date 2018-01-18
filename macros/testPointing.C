double testPointing(double dtheta = 0, int run = 332, int event = 55448680)
// double testPointing(double dtheta = 0, int run = 352, int event = 60832108)   
{

  AnitaDataset d(run); 
  d.getEvent(event); 

  UsefulAdu5Pat gps(d.gps()); 

  double wais_phi, wais_theta; 
  gps.getThetaAndPhiWaveWaisDivide(wais_theta, wais_phi); 
  wais_theta += dtheta * TMath::DegToRad(); 

  double rampdem_wais_alt = AnitaLocations::getWaisAltitude(); 
  printf("Real WAIS lat/lon/alt: %g %g %g\n", AnitaLocations::getWaisLatitude(), AnitaLocations::getWaisLongitude(), AnitaLocations::getWaisAltitude()); 
  printf(" According to AnitaGeomTool, WAIS at phi=%g theta=%g\n", wais_phi* TMath::RadToDeg(), wais_theta* TMath::RadToDeg()); 

  PayloadParameters pp0(d.gps(), AntarcticCoord(AntarcticCoord::WGS84, AnitaLocations::getWaisLatitude(), AnitaLocations::getWaisLongitude(), rampdem_wais_alt )); 
  printf(" According to PayloadParameters, wais at phi=%g, theta=%g\n", pp0.source_phi, pp0.source_theta); 


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


  PayloadParameters pp3; 
  int success = PayloadParameters::findSourceOnContinent(pp0.source_theta, pp0.source_phi, d.gps(), &pp3); 

  if (success==1) 
  {
    AntarcticCoord c3 = pp3.source.as(AntarcticCoord::WGS84); 
    printf("findSourceOnContinent projects to lat/lon/alt: %g %g %g\n",c3.x,c3.y,c3.z); 
    printf("Payload parameters thinks it's at phi=%g, theta=%g\n", pp3.source_phi, pp3.source_theta); 
  }
  else
  {
    printf("findSourceOnContinent found no solution :(\n"); 

  }

  double x0 = pp0.payload.as(AntarcticCoord::STEREOGRAPHIC).x; 
  double y0 = pp0.payload.as(AntarcticCoord::STEREOGRAPHIC).y; 


  TFile * fout = new TFile("refraction.root","RECREATE"); 
  TH2I * old_h = new TH2I("traceBack","TraceBack", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3);
  TH2I * new_h3 = new TH2I("traceBack3","TraceBack3", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3);   
  TH2I * new_h = new TH2I("findSource","FindSource", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3); 
  TH2I * new_h_col = new TH2I("findSourceColl","FindSourceColl", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3); 
  TH2I * new_h_ref = new TH2I("findSourceRef","FindSourceRefract", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3); 
  TH2 * old_payload_el = new TProfile2D("traceBackPayloadEl","TraceBackPayloadEl", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3);
  TH2 * new_h3_payload_el = new TProfile2D("traceBack3PayloadEl","TraceBack3PayloadEl", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3);   
  TH2 * new_payload_el = new TProfile2D("findSourcePayloadEl","findSourcePayloadEl", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3); 
  TH2 * new_payload_el_col = new TProfile2D("findSourcePayloadElCol","findSourcePayloadElCol", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3); 
  TH2 * new_payload_el_ref = new TProfile2D("findSourcePayloadElRef","findSourcePayloadElRef", 1000, x0 - 800e3, x0+800e3, 1000, y0 - 800e3, y0+800e3); 

  Refraction::SphRay m; 
  for (double phi = 0; phi <=360; phi+= 0.2) 
  {
    printf("%g\n",phi); 
    double step = 0.001;
    for (double sin_theta = 0.0; sin_theta < 0.9; sin_theta+= step)
    {
        double theta = asin(sin_theta) * TMath::RadToDeg(); 
        trace = gps.traceBackToContinent(phi*TMath::DegToRad(), theta*TMath::DegToRad(), &lon, &lat,&alt,0); 
        if (trace) 
        {
          AntarcticCoord w(AntarcticCoord::WGS84,lat,lon,alt); 
          w.to(AntarcticCoord::STEREOGRAPHIC); 
          old_h->Fill(w.x, w.y); 
          PayloadParameters pw(d.gps(), w); 
          old_payload_el->Fill(w.x,w.y, pw.payload_el); 

        }


        trace = gps.traceBackToContinent3(phi*TMath::DegToRad(), theta*TMath::DegToRad(), &lon, &lat,&alt,0); 
        if (trace)
        {
          AntarcticCoord w(AntarcticCoord::WGS84,lat,lon,alt);
          w.to(AntarcticCoord::STEREOGRAPHIC);
          new_h3->Fill(w.x, w.y);
          PayloadParameters pw(d.gps(), w);
          new_h3_payload_el->Fill(w.x,w.y, pw.payload_el);
        }

        PayloadParameters::findSourceOnContinent(theta,phi, d.gps(), &pp3); 
        AntarcticCoord w = pp3.source.as(AntarcticCoord::STEREOGRAPHIC); 
        new_h->Fill(w.x,w.y); 
        new_payload_el->Fill(w.x,w.y, pp3.payload_el); 

        PayloadParameters::findSourceOnContinent(theta,phi, d.gps(), &pp3,0, 50); 
        w = pp3.source.as(AntarcticCoord::STEREOGRAPHIC); 
        new_h_col->Fill(w.x,w.y); 
        new_payload_el_col->Fill(w.x,w.y, pp3.payload_el); 

        PayloadParameters::findSourceOnContinent(theta,phi, d.gps(), &pp3,&m); 
        w = pp3.source.as(AntarcticCoord::STEREOGRAPHIC); 
        new_h_ref->Fill(w.x,w.y); 
        new_payload_el_ref->Fill(w.x,w.y, pp3.payload_el); 


    }
  }

  TCanvas * cc = new TCanvas; 
  cc->Divide(5,2);
  cc->cd(1);
  old_h->Draw("colz"); 
  cc->cd(2);
  new_h->Draw("colz"); 
  cc->cd(3);
  new_h_col->Draw("colz"); 
  cc->cd(4);
  new_h_ref->Draw("colz"); 
  cc->cd(5);
  new_h3->Draw("colz");

  cc->cd(6);
  old_payload_el->Draw("colz");
  cc->cd(7);
  new_payload_el->Draw("colz");
  cc->cd(8);
  new_payload_el_col->Draw("colz"); 
  cc->cd(9);
  new_payload_el_ref->Draw("colz");
  cc->cd(10);
  new_h3_payload_el->Draw("colz");


  cc->SaveAs("refraction.pdf"); 

  fout->cd(); 
  old_h->Write(); 
  new_h->Write(); 
  new_h_col->Write(); 
  new_h_ref->Write(); 
  old_payload_el->Write(); 
  new_payload_el->Write(); 
  new_payload_el_col->Write(); 
  new_payload_el_ref->Write(); 


  return wais_phi * TMath::RadToDeg() - pp0.source_phi; 

}



