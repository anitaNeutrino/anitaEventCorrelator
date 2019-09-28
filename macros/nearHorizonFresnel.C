

#define USE_SPLINE

const double horizon = -5.93934873; 
const double above = 0.02; 
const double delta = 1e-2; 

void nearHorizonFresnel(double start_alt = 38.576e3, double surface_height=2580, 
    double direct_angle = horizon+above, double dthrown = above*above*delta, double first_extra_skip = above-2*delta*above*above, 
    int n_below_horizon=50)
{


  AntarcticAtmosphere::ExponentialRefractivity  m(313, 1./8000);
  Refraction::RaytracerSpherical::Setup s;
  Refraction::RaytracerSpherical::Result r; 
  s.start_alt = start_alt; 
  s.end_alt = 80e3; 
  s.surface_alt = surface_height; 

  Refraction::RaytracerSpherical ray(&m); 
  ray.step_size=10; 


  std::vector<double> thrown;
  std::vector<double> reflect_angle;
  std::vector<TGraph *> xy;
  std::vector<TGraph *> xt;
  std::vector<double> X; 
  std::vector<double> S; 
  std::vector<bool> drawme; 

  drawme.push_back(true); 

  double th = direct_angle; 
  while (thrown.size() <=n_below_horizon)
  {
    s.thrown_payload_angle = th;
    printf("Throwing %f\n",th);
    ray.raytrace(&s,&r); 
    if (th == direct_angle) 
    { 
      if (r.reflected)
      {
        fprintf(stderr,"%g isn't a direct ray!\n", th); 
        return; 
      }

      X.resize(ray.nSteps()); 
      S.resize(ray.nSteps()); 
      printf("n steps: %d\n", ray.nSteps()); 
      memcpy(&X[0], ray.lastGrammage(), sizeof(double) * X.size()); 
      memcpy(&S[0], ray.lastS(), sizeof(double) * S.size()); 
    }


    if (th == direct_angle || r.reflected)
    {
      xy.push_back(ray.makeXYGraph());
      xt.push_back(ray.makeXTGraph());
      thrown.push_back(th); 
      reflect_angle.push_back(r.reflection_angle); 
    }
    if (r.reflected)
    {
      printf("...reflected! angle: %g\n", r.reflection_angle); 
    }
    if (th == direct_angle) th-=first_extra_skip; 
    th-=dthrown; 
  }

  //now let's spline the thrown xy and xt; 

 
#ifdef USE_SPLINE
  TSpline3 sxy0("",xy[0]->GetX(), xy[0]->GetY(), xy[0]->GetN()); 
  TSpline3 sxt0("",xt[0]->GetX(), xt[0]->GetY(), xt[0]->GetN()); 
#else
  TGraph & sxy0 = *xy[0]; 
  TGraph & sxt0 = *xt[0]; 
  sxy0.SetBit(TGraph::kIsSortedX) ;
  sxt0.SetBit(TGraph::kIsSortedX) ;
#endif

  TGraph * direct_grammage = new TGraph(xy[0]->GetN());
  direct_grammage->SetTitle(Form("grammage (%g)", above)); 
  direct_grammage->GetXaxis()->SetTitle("x (m)"); 
  direct_grammage->GetYaxis()->SetTitle("X (g/cm^{2})"); 
  TGraph *direct_grammage_alt = new TGraph;
  direct_grammage_alt->SetTitle(Form("Interaction altitude vs. grammage (%g)", above)); 
  direct_grammage_alt->GetYaxis()->SetTitle("interaction alt (m)"); 
  direct_grammage_alt->GetXaxis()->SetTitle("X (g/cm^{2})"); 
  TGraph *direct_grammage_dist = new TGraph;
  direct_grammage_dist->SetTitle(Form("path distance vs. grammage (%g)", above)); 
  direct_grammage_dist->GetYaxis()->SetTitle("path distance (m)"); 
  direct_grammage_dist->GetXaxis()->SetTitle("X (g/cm^{2})"); 



  for (int i = 0; i < xy[0]->GetN(); i++) 
  {

    direct_grammage->SetPoint(i, xy[0]->GetX()[i], X[i]); 
    direct_grammage_dist->SetPoint(i, X[i], S[i]); 

    double R2 = pow(xy[0]->GetX()[i],2) + pow(xy[0]->GetY()[i],2);
    if (i > 1)
    {
      double last_R2 = pow(xy[0]->GetX()[i-1],2) + pow(xy[0]->GetY()[i-1],2);
      if (R2 > last_R2) 
      {
        direct_grammage_alt->SetPoint(direct_grammage_alt->GetN(), X[i],sqrt(R2)-s.R_c);
      }
    }

  }

  direct_grammage->SetBit(TGraph::kIsSortedX) ;;



  TGraph *crossing_x = new TGraph; 
  crossing_x->SetTitle(Form("Crossing x (%g)", above)); 
  crossing_x->GetXaxis()->SetTitle("Reflected ray angle (deg)"); 
  crossing_x->GetYaxis()->SetTitle("x (m)"); 

  TGraph *crossing_alt = new TGraph; 

  crossing_alt->SetTitle(Form("Crossing alt (%g)", above)); 
  crossing_alt->GetXaxis()->SetTitle("Reflected ray angle (deg)"); 
  crossing_alt->GetYaxis()->SetTitle("alt (m)"); 


  TGraph *crossing_dt = new TGraph; 

  crossing_dt->SetTitle(Form("Crossing dt (%g)", above)); 
  crossing_dt->GetXaxis()->SetTitle("Reflected ray angle (deg)"); 
  crossing_dt->GetYaxis()->SetTitle("dt (ns)"); 

  TGraph *crossing_dtheta = new TGraph; 
  crossing_dtheta->SetTitle(Form("Crossing d#theta (%g)", above)); 
  crossing_dtheta->GetXaxis()->SetTitle("Reflected ray angle (deg)"); 
  crossing_dtheta->GetYaxis()->SetTitle("d#theta (deg)"); 



  TGraph *grammage_dt = new TGraph; 

  grammage_dt->SetTitle(Form("dt vs grammage (%g)", above)); 
  grammage_dt->GetXaxis()->SetTitle("X (g/cm^2)"); 
  grammage_dt->GetYaxis()->SetTitle("dt (ns)"); 

  TGraph * grammage_dtheta = new TGraph;
  grammage_dtheta->SetTitle(Form("d#theta_{obs} vs grammage (%g)", above)); 
  grammage_dtheta->GetXaxis()->SetTitle("X (g/cm^2)"); 
  grammage_dtheta->GetYaxis()->SetTitle("d#theta_{obs} (deg)"); 



  for (int i = 1; i < thrown.size(); i++) 
  {
    bool dodraw = false;
   //check if there's a crossing or not 
    for (int j = 2; j < xy[i]->GetN(); j++) 
    {

      if ( (xy[i]->GetY()[j] < sxy0.Eval(xy[i]->GetX()[j])) != (xy[i]->GetY()[j-1] < sxy0.Eval(xy[i]->GetX()[j-1])))
      {
        //we found a crossing! now let's find the root 

#ifdef USE_SPLINE
        TSpline3 sxy("",xy[i]->GetX(), xy[i]->GetY(), xy[i]->GetN()); 
        TSpline3 sxt("",xt[i]->GetX(), xt[i]->GetY(), xt[i]->GetN()); 
#else
        TGraph & sxy = *xy[i]; 
        TGraph & sxt = *xt[i]; 
        sxy.SetBit(TGraph::kIsSortedX) ;;
        sxt.SetBit(TGraph::kIsSortedX) ;;
#endif

        TF1 root_fn("root_fn", [&](double * x, double * p) {  return sxy.Eval(*x)-sxy0.Eval(*x); }, xy[i]->GetX()[j-1], xy[i]->GetX()[j], 0); 
        double root = root_fn.GetX(0); 

        crossing_x->SetPoint(crossing_x->GetN(), thrown[i], root);
        double y = sxy.Eval(root); 
        double alt = sqrt(root*root+y*y)-s.R_c;
        printf("Found crossing for reflected ray %g at %g,%g (alt=%g)\n", thrown[i], root,y,alt); 
        dodraw=true; 
        crossing_alt->SetPoint(crossing_alt->GetN(), thrown[i], alt);
        double dt = 1e9*(sxt0.Eval(root) - sxt.Eval(root));
        crossing_dt->SetPoint(crossing_dt->GetN(), thrown[i], dt);

        //TODO, throw another ray here based on off-cone angles? 
        double grammage = direct_grammage->Eval(root); 
        grammage_dt->SetPoint(grammage_dt->GetN(), grammage, dt); 

        //estimate the delta_t , we need to find the right segment of the original one
        TVector3 d3(1, sxy0.Derivative(root),0); 
        TVector3 r3(1, sxy.Derivative(root),0); 

        double dtheta = TMath::RadToDeg()*d3.Angle(r3); 
        crossing_dtheta->SetPoint(crossing_dtheta->GetN(), thrown[i], dtheta);
        grammage_dtheta->SetPoint(grammage_dtheta->GetN(), grammage, dtheta);

      }
    }
    drawme.push_back(dodraw); 
  }

  TCanvas * cx = new TCanvas("cx","x",1000,600);
  cx->Divide(2,1); 
  cx->cd(1); 
  direct_grammage->Draw("alp");
  cx->cd(2); 
  crossing_x->Draw("alp");

  TCanvas * calt = new TCanvas("calt","alt",1000,600);
  calt->Divide(2,1); 
  calt->cd(1); 
  direct_grammage_alt->Draw("alp");
  calt->cd(2); 
  crossing_alt->Draw("alp");


  TCanvas * cdt = new TCanvas("cdt","dt",1000,600);
  cdt->Divide(2,1); 
  cdt->cd(1); 
  grammage_dt->Draw("alp");
  cdt->cd(2); 
  crossing_dt->Draw("alp");

  TCanvas * cdtheta = new TCanvas("cdtheta","dtheta",1000,600);
  cdtheta->Divide(2,1); 
  cdtheta->cd(1); 
  grammage_dtheta->Draw("alp");
  cdtheta->cd(2); 
  crossing_dtheta->Draw("alp");


  /*

  TCanvas * cdraw = new TCanvas("cdraw","draw",1000,600); 

  TMultiGraph * mg = new TMultiGraph;
  for (int i = 0; i < thrown.size(); i++) 
  {
    xy[i]->SetTitle(Form("%g deg", thrown[i]));
    if (i > 0 && drawme[i])
    {
      mg->Add(xy[i]); 
    }
  }

  mg->Draw("a plc pmc"); 

  TGraph * earth = new TGraph(900); 
  TGraph * galts[5] = {0};
  double alts[5] = { 5000,7500,10000,12500,15000}; 
  double Rs[5];
  for (int i = 0; i < 5; i++) 
  {
    Rs[i] = s.R_c + alts[i]; 
    galts[i] = new TGraph(900); 
  }

  double R = s.R_c + surface_height;
  for (int i = 0; i < 900; i++) 
  {
    double ang = i/10. * TMath::DegToRad();
    earth->SetPoint(i, R*cos(ang), R*sin(ang)); 
    for (int j = 0; j < 5; j++) 
    {
      galts[j]->SetPoint(i,Rs[j]*cos(ang), Rs[j]*sin(ang)); 
    }
  }
  earth->SetTitle("Earth Surface"); 
  earth->SetLineColor(1); 
  earth->Draw("lsame"); 

  for (int i = 0; i < 5; i++) 
  {
    galts[i]->SetLineStyle(2+i); 
    galts[i]->Draw("lsame"); 
    galts[i]->SetTitle(Form("alt = %g m",alts[i])); 
  }

  cdraw->BuildLegend(); 
  */

  TCanvas * cfresnel = new TCanvas("cfresnel","fresnel",1800,1000); 

  cfresnel->SetGridx();
  cfresnel->SetGridy();
  TGraph * fresnel = new TGraph; 
  fresnel->SetTitle(Form("Fresnel loss for #theta=%g (X=800 g/cm^{2})",above)); 

  fresnel->GetXaxis()->SetTitle("Frequency (GHz)"); 
  fresnel->GetYaxis()->SetTitle("Fresnel Loss or Gain (dB)"); 
  double lag = grammage_dt->Eval(800); 
  for (double f = 0.2; f < 0.8; f+=0.001) 
  {
    double lambda = 1. / f; 
    fresnel->SetPoint(fresnel->GetN(), f,  20*log10(  1 - cos(2*TMath::Pi()*fabs(lag)/lambda)));
  }

  fresnel->SetLineWidth(2); 
  fresnel->Draw(); 





  TCanvas *csave = new TCanvas("csave","save",1800,1000); 
  csave->Divide(2,2); 
  csave->cd(1); 
  gPad->SetGridx();
  gPad->SetGridy();
  grammage_dt->SetLineWidth(2); 
  grammage_dt->GetXaxis()->SetRangeUser(0,2000); 
  grammage_dt->Draw("alp");
  csave->cd(2); 
  gPad->SetGridx();
  gPad->SetGridy();
  grammage_dtheta->SetLineWidth(2); 
  grammage_dtheta->GetXaxis()->SetRangeUser(0,2000); 
  grammage_dtheta->Draw("alp");
  csave->cd(3); 
  gPad->SetGridx();
  gPad->SetGridy();
  direct_grammage_alt->SetLineWidth(2); 
  direct_grammage_alt->GetXaxis()->SetRangeUser(0,2000); 
  direct_grammage_alt->Draw("alp");
  csave->cd(4); 
  gPad->SetGridx();
  gPad->SetGridy();
  direct_grammage_dist->GetXaxis()->SetRangeUser(0,2000); 
  direct_grammage_dist->SetLineWidth(2); 
  direct_grammage_dist->Draw("alp");







}
