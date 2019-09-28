


double the_max_phi = 16; 
AntarcticAtmosphere::AtmosphericModel * makeAtmosphere(int nsoundings, double max_phi=16, std::vector<AntarcticAtmosphere::AtmosphericModel *> * store = 0,
     int * years = 0, int * months = 0, int * days = 0, int * earlies = 0) 
{


  AntarcticAtmosphere::InterpolatedAtmosphere * m = new AntarcticAtmosphere::InterpolatedAtmosphere; 

  double dphi = max_phi/(nsoundings-1); 
  for (int i = 0; i < nsoundings; i++) 
  {


    int year = 2015+gRandom->Integer(4); 
    int month = 12; 
    int day = 1 + gRandom->Integer(31); 
    bool early = gRandom->Integer(2); 

    if (days) days[i] = day; 
    if (years) years[i] = year; 
    if (months) months[i] = month; 
    if (earlies) earlies[i] = early; 

    AntarcticAtmosphere::SPRadiosonde * r = new AntarcticAtmosphere::SPRadiosonde(year,month,day,early); 
    if (store) store->push_back(r); 
    m->addModel(r, nsoundings > 1 ? dphi*i : 0); 
  }

  return m; 
}


int find_crossings(Refraction::RaytracerSpherical * r, TSpline3 * sxy, std::vector<int> * xs)
{

  int ncross = 0;
  for (int j =r->nSteps()-1; j > 1; j--)
  {
    if (r->lastX()[j] > sxy->GetXmax()) continue; 
    if (r->lastY()[j]  < sxy->Eval(r->lastX()[j]) != (r->lastY()[j-1] < sxy->Eval(r->lastX()[j-1])))
    {
      xs->push_back(j-1); 
      ncross++;
    }
  }

  return ncross; 
}







void radiosondeRaytrace(double direct_angle = -5.9,  double dtheta=1e-5,  int nsoundings=5, int natm = 50, double start_alt = 39e3, double surface_height=2900) 
{



  //figure out an output file name that's not in use
  //
  int ifile = 0; 
  TFile * f = 0;
  while (!f || f->IsZombie() || !f->IsOpen()) 
  {
    f = new TFile(Form("radiosonde_out/radiosonde_raytrace_%d_%g.%d.root", nsoundings, direct_angle,ifile++),"CREATE"); 
  }

  gRandom->SetSeed(1+ifile); 

  TTree * t = new TTree("intersection_tree","Radiosonde raytracing"); 
  TTree * tdirect = new TTree("direct_tree","Radiosonde raytracing direct tree"); 
  TTree * tplot= new TTree("plot_tree","Radiosonde raytracing"); 

  TH2 * hN = new TH2F("hN","Refractivity;phi;alt", 100, 0, the_max_phi, (50e3-surface_height)/5,surface_height,50e3); 
  TGraph * gref = new TGraph; 
  TMultiGraph * gdirect = new TMultiGraph; 
  TMultiGraph * greflect = new TMultiGraph; 

  tplot->Branch("hN",&hN);
  tplot->Branch("ref",&gref);
  tplot->Branch("gdirect",&gdirect);
  tplot->Branch("greflect",&greflect);


  t->Branch("direct_angle",&direct_angle); 
  tdirect->Branch("direct_angle",&direct_angle); 

  double intersection_X; 
  double intersection_dt; 
  double intersection_S; 
  double intersection_dS; 
  double intersection_dtheta; 
  double intersection_alt; 
  double thrown;
  int reflected;
  int ncrossings = 0; 
  int icross = 0; 

  int sounding_day[nsoundings];
  int sounding_year[nsoundings];
  int sounding_month[nsoundings];
  int sounding_early[nsoundings];

  int iatm; 
  t->Branch("intersection_X", &intersection_X); 
  t->Branch("iatm", &iatm); 
  tdirect->Branch("iatm", &iatm); 
  t->Branch("intersection_S", &intersection_S); 
  t->Branch("intersection_dS", &intersection_dS); 
  t->Branch("intersection_dt", &intersection_dt); 
  t->Branch("intersection_dtheta", &intersection_dtheta); 
  t->Branch("intersection_alt", &intersection_alt); 
  t->Branch("thrown_angle", &thrown); 
  t->Branch("ncrossings", &ncrossings); 
  t->Branch("icross", &icross); 
  t->Branch("reflected", &reflected); 
  tdirect->Branch("reflected", &reflected); 
  t->Branch("nsoundings", &nsoundings); 
  t->Branch("sounding_day[nsoundings]", &sounding_day,"sounding_day[nsoundings]/I"); 
  t->Branch("sounding_month[nsoundings]", &sounding_month,"sounding_month[nsoundings]/I"); 
  t->Branch("sounding_year[nsoundings]", &sounding_year,"sounding_year[nsoundings]/I"); 
  t->Branch("sounding_early[nsoundings]", &sounding_early,"sounding_early[nsoundings]/I"); 
  tdirect->Branch("nsoundings", &nsoundings); 
  tdirect->Branch("sounding_day[nsoundings]", &sounding_day,"sounding_day[nsoundings]/I"); 
  tdirect->Branch("sounding_month[nsoundings]", &sounding_month,"sounding_month[nsoundings]/I"); 
  tdirect->Branch("sounding_year[nsoundings]", &sounding_year,"sounding_year[nsoundings]/I"); 
  tdirect->Branch("sounding_early[nsoundings]", &sounding_early,"sounding_early[nsoundings]/I"); 
  double alt_min; 
  t->Branch("alt_min",&alt_min);
  tdirect->Branch("alt_min",&alt_min);



  Refraction::RaytracerSpherical::Setup s;
  s.start_alt = start_alt; 
  s.end_alt = 80e3; 
  s.surface_alt = surface_height; 
  s.thrown_payload_angle = direct_angle; 
  Refraction::RaytracerSpherical::Result r; 
  for (iatm = 0; iatm < natm; iatm++)
  {
    //randomly generate the right number of atmospheres... 
    std::vector<AntarcticAtmosphere::AtmosphericModel*> store; 
    AntarcticAtmosphere::AtmosphericModel * m = makeAtmosphere(nsoundings,the_max_phi,&store, sounding_year, sounding_month, sounding_day, sounding_early);

    printf("Filling in atmosphere..\n"); 
    for (int i = 1; i <=hN->GetNbinsX(); i++) 
    {

     double phi= hN->GetXaxis()->GetBinCenter(i);
      for (int j = 1; j <=hN->GetNbinsY(); j++) 
      {
        double alt= hN->GetYaxis()->GetBinCenter(j);
        hN->SetBinContent(i,j, m->get(alt, AntarcticAtmosphere::REFRACTIVITY, phi)); 
      }

    }
    printf("..done\n");

    Refraction::RaytracerSpherical ray(m); 
    s.thrown_payload_angle = direct_angle; 
    //first lets do the direct path 
    ray.raytrace(&s,&r); 

    reflected = r.reflected; 

    if (reflected) 
    {
      printf("reference angle %g reflected on this trial!", direct_angle); 
      alt_min = surface_height; 
      ray.makePhiAltGraph(gref);
      tplot->Fill();
      tdirect->Fill(); 
      continue; 
    }

    //otherwise, let's find the minimum altitude for the  reference path
    alt_min = ray.lastMinAlt(); 
    tdirect->Fill(); 
    
    //store reference path 
    TGraph grammage(ray.nSteps(), ray.lastX(), ray.lastGrammage()); 
    grammage.SetBit(TGraph::kIsSortedX) ;


    TSpline3 sxy_ref("",(double*) ray.lastX(), (double*) ray.lastY(), ray.nSteps());
    TSpline3 sxt_ref("",(double*) ray.lastX(), (double*) ray.lastT(), ray.nSteps());
    TSpline3 sxs_ref("",(double*) ray.lastX(), (double*) ray.lastS(), ray.nSteps());

    //store the paty

    ray.makePhiAltGraph(gref);
    //raytrace a bit less
    s.end_alt = 40e3; 

    bool going_up = true; 
    double real_dtheta =dtheta; 
    while(true) 
    {
      s.thrown_payload_angle += going_up ? dtheta : -dtheta; 
      printf("Throwing %g\n", s.thrown_payload_angle); 
      ray.raytrace(&s,&r); 

      //now let's look for a crossing

      std::vector<int> crossings; 
      ncrossings = find_crossings(&ray, &sxy_ref, &crossings); 

      //no crossing 
      if (!ncrossings) 
      {
        if (going_up) 
        {
          printf("Going down at %g\n", s.thrown_payload_angle); 
          going_up = false; 
          s.thrown_payload_angle = direct_angle; 
        }
        continue; 
      }
      //otherwise, let's find where we cross

      bool should_save = false; 

      //loop over all crossings 

      std::vector<double> Xs; 
      std::vector<double> dts; 
      thrown = s.thrown_payload_angle; 
      reflected = r.reflected; 
      for (icross= 0; icross< ncrossings; icross++) 
      {
        int cross_index = crossings[icross]; 

        int min_i = TMath::Max(0,cross_index-10); 
        int max_i = TMath::Min(ray.nSteps(),cross_index+10); 

        TSpline3 sxy("",(double*) ray.lastX()+min_i, (double*) ray.lastY()+min_i, max_i-min_i+1);
        TSpline3 sxt("",(double*) ray.lastX()+min_i, (double*) ray.lastT()+min_i, max_i-min_i+1);
        TSpline3 sxs("",(double*) ray.lastX()+min_i, (double*) ray.lastS()+min_i, max_i-min_i+1);

        TF1 root_fn("root_fn", [&](double * x, double * p) {  return sxy.Eval(*x)-sxy_ref.Eval(*x); }, ray.lastX()[cross_index], ray.lastX()[cross_index+1], 0); 
        double root = root_fn.GetX(0); 

        intersection_X = grammage.Eval(root); 
        if (intersection_X > 600 && intersection_X < 1000) should_save = true; 
        intersection_S = sxs_ref.Eval(root); 
        intersection_dS =sxs.Eval(root) -  intersection_S;
        intersection_dt = 1e9*(sxt.Eval(root) - sxt_ref.Eval(root)); 
        dts.push_back(intersection_dt);
        Xs.push_back(intersection_X); 
        printf("Crossing! X= %g, dt=%g, reflected=%d\n", intersection_X, intersection_dt, reflected); 

        TVector3 d3(1, sxy_ref.Derivative(root),0); 
        TVector3 r3(1, sxy.Derivative(root),0); 
        intersection_dtheta = TMath::RadToDeg()*d3.Angle(r3); 
        double y = sxy.Eval(root); 
        intersection_alt = sqrt(root*root + y*y) - s.R_c; 
        alt_min = ray.lastMinAlt(); 

        t->Fill(); 
      }

      if (should_save) 
      {
        TGraph * gpath = ray.makePhiAltGraph(); 
        TString title = Form("thrown=%g%s", thrown, reflected?" (reflected) ":" ");

        if (Xs.size()==1) title += Form("X=%g", Xs[0]);
        else 
        {
          title += Form("Xs=%g,",Xs[0]); 
          for (int iX = 1;iX<Xs.size();iX++) title+=Form(",%g",Xs[iX]);
        }

        if (dts.size()==1) title += Form("dt=%g", dts[0]);
        else 
        {
          title += Form("dts=%g,",dts[0]); 
          for (int idt = 1;idt<dts.size();idt++) title+=Form(",%g",dts[idt]);
        }

 

        gpath->SetTitle(title);
        if (reflected) greflect->Add(gpath);
        else gdirect->Add(gpath);
      }


      //break when reflected rays have unrealistic X or we already have 50 reflected graphs... 
     
      if (!going_up &&  reflected && intersection_X > 1e3) 
      {
        //first make sure we actually have some reflected graphs.. 
        if (!greflect->GetListOfGraphs())
        {
          printf("Adjusting to try to get more near horizon\n"); 
          s.thrown_payload_angle += dtheta; 
          dtheta*=0.1; 
        }
        else
        {
          break; 
        }

      }


      //if we already have more than 20 reflected graphs in our range, we can make our step bigger
      if (greflect->GetListOfGraphs() && greflect->GetListOfGraphs()->GetEntries() > 20) dtheta *=2; 

    }

    dtheta = real_dtheta; 
    printf("Stopped at %g\n", s.thrown_payload_angle); 

    tplot->Fill(); 


    delete gdirect; 
    delete greflect; 
    gdirect = new TMultiGraph; 
    greflect = new TMultiGraph; 


    for (auto im : store) delete im; 
    delete m; 
  }


  t->Write(); 
  tdirect->Write();
  tplot->Write(); 

  delete f; 
}
 

