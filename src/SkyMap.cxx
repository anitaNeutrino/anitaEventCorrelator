#include "SkyMap.h" 
#include "TMath.h" 
#include "TPaletteAxis.h" 
#include "TEllipse.h" 
#include "FFTtools.h" 
#include "TPad.h" 


static void convert(TGraph * g) 
{
  for (int i = 0; i < g->GetN(); i++)
  {
    double x,y; 
    SkyMap::toMollweide(g->GetX()[i],g->GetY()[i], x,y); 
    g->GetX()[i] = x; 
    g->GetY()[i] = y; 
  }

}

static void paintSpecials(bool reverse = true)
{

  static std::vector<TObject*> curves; 

  if (curves.size() == 0) 
  {
    //Make the ellipse
    TEllipse * ell =  new TEllipse(0,0,2,1); 
    curves.push_back(ell); 
    ell->SetLineColor(16); 
    ell->SetLineWidth(3); 
    ell->SetFillStyle(0); 


    // latitudes 
    for (int ilat = -60; ilat <= 60; ilat += 30)
    {
      TGraph * glat = new TGraph(3); 
      glat->SetPoint(0,-179.9999,ilat); 
      glat->SetPoint(1,0,ilat); 
      glat->SetPoint(2,179.9999,ilat); 
      glat->SetLineColor(16); 
      glat->SetLineWidth(2); 
      convert(glat); 
      curves.push_back(glat); 
    }

    //longitudes

    for (int ilong = -150; ilong <=150; ilong+=30)
    {
      TGraph * glon = new TGraph(181); 
      for (int i = 0; i <= 180; i++) glon->SetPoint(i, ilong, i-90); 
      convert(glon); 
      glon->SetLineColor(16); 
      glon->SetLineWidth(2); 
      curves.push_back(glon); 
    }
  }


  for (unsigned i = 0; i < curves.size(); i++) curves[i]->Paint(); 

}




static int skymap_counter = 0; 
SkyMap::SkyMap(double lon_0, int nbinsx , int nbinsy , const TH2 * background, const std::vector<const TMarker * > * markers, const std::vector<const TGraph *> * graphs )
  : reverseX(true), lon_0(lon_0), sky_background(Form("skymap_hist_%d", skymap_counter++),  "Sky Map", nbinsx, -2, 2., nbinsy*1.1, -1, 1), 
    left(-2.1,0, TString::Format("%g",reverseX ? lon_0 + 180: lon_0-180)), 
    right(2.1,0, TString::Format("%g",reverseX?lon_0 -180:lon_0+180)) 
{
  zaxis = 0; 
  range_set = false; 
  sky_background.SetDirectory(0); 
  left.SetTextAlign(32); 
  right.SetTextAlign(12); 
  left.SetTextColor(13); 
  right.SetTextColor(13); 

  if(background) setBackground(background); 
  sky_background.SetStats(false); 

  if (markers) 
  {
    for (unsigned i = 0; i < markers->size(); i++) addMarker(markers->at(i)); 
  }

  if (graphs) 
  {
    for (unsigned i = 0; i < graphs->size(); i++) addGraph(graphs->at(i)); 
  }
}


static double aux_angle(double lat_rad, double eps = 1e-12, int maxiter = 1000) 
{
  if (cos(lat_rad) < 0.0001) return lat_rad; 

  double two_theta = lat_rad; 
  int niter = 0; 
  while(true) 
  {
    double d_two_theta = -( two_theta + sin (two_theta) - TMath::Pi() * sin(lat_rad)) / (1 + cos(two_theta));  

    two_theta += d_two_theta; 
    if (niter > 5 && (fabs(d_two_theta) < eps || niter >= maxiter)) return two_theta/2; 
    niter++; 
  }
}

void SkyMap::toMollweide(double lon, double lat, double &x, double &y, double lon_0, bool revx)
{
  double lon_rad = FFTtools::wrap(lon-lon_0,360,0) * TMath::DegToRad(); 
  double lat_rad = lat * TMath::DegToRad(); 
  double theta = aux_angle(lat_rad); 
  x = 2/ TMath::Pi() * lon_rad * cos(theta); 
  if (revx) x =-x; 
  y =  sin(theta); 
}

void SkyMap::fromMollweide(double x, double y, double &lon, double & lat, double lon_0, bool revx) 
{
  if (revx) x = -x; 
  double theta = asin(y); 
  double lat_rad = asin(  (2*theta+sin(2*theta))/TMath::Pi()  ); 
  double lon_rad = TMath::Pi() * x / (2 * cos(theta)); 
  lat = TMath::RadToDeg()*lat_rad; 
  lon = FFTtools::wrap(TMath::RadToDeg()*lon_rad + lon_0,360,0); 
}


void SkyMap::setBackground(const TH2 * bg, bool ra_in_hours) 
{
  if (!bg) 
  {
    sky_background.Reset(); 
    if (zaxis) delete zaxis; 
    zaxis = 0; 
  }

  double scale = 1; 

  double lat,lon; 
  for (int i = 1; i <= sky_background.GetNbinsX(); i++) 
  {
    double x = sky_background.GetXaxis()->GetBinCenter(i); 
    for (int j = 1; j <=sky_background.GetNbinsY(); j++) 
    {
      double y = sky_background.GetYaxis()->GetBinCenter(j); 
      if ((x*x)/4 + y*y >1) continue; 
      fromMollweide(x,y,lon,lat,lon_0, reverseX); 
      if (ra_in_hours) lon/=15; 

      //convert latitude/lon to histogram units 
      lon = FFTtools::wrap(lon, bg->GetXaxis()->GetXmax() - bg->GetXaxis()->GetXmin(), 0.5 * (bg->GetXaxis()->GetXmin() + bg->GetXaxis()->GetXmax())); 
      lat = FFTtools::wrap(lat, bg->GetYaxis()->GetXmax() - bg->GetYaxis()->GetXmin(), 0.5 * (bg->GetYaxis()->GetXmin() + bg->GetYaxis()->GetXmax())); 
      sky_background.SetBinContent(i,j, ((TH2*) bg)->Interpolate(lon,lat)); 
    }
  }
  zaxis =  new TPaletteAxis(2.0,0.2,2.05,0.9, &sky_background); 
}

void SkyMap::addMarker(const TMarker * marker) 
{
  double x,y; 
  toMollweide(marker->GetX(), marker->GetY(),x,y, lon_0, reverseX);; 
  sky_markers.push_back(TMarker(*marker)); 
  sky_markers[sky_markers.size()-1].SetX(x); 
  sky_markers[sky_markers.size()-1].SetY(y); 
}

void SkyMap::addText(const TLatex * text) 
{
  double x,y; 
  toMollweide(text->GetX(), text->GetY(),x,y, lon_0, reverseX);; 
  sky_text.push_back(TLatex(*text)); 
  sky_text.back().SetX(x); 
  sky_text.back().SetY(y); 
}

void SkyMap::addGraph(const TGraph * gr, bool ra_in_hours) 
{
  sky_graphs.push_back(TGraph(gr->GetN())); 
  TGraph & g = sky_graphs[sky_graphs.size()-1]; 
  double xscale = 1; 
  if (ra_in_hours) xscale = 15; 

  g.SetLineStyle(gr->GetLineStyle()); 
  g.SetLineColor(gr->GetLineColor()); 
  g.SetLineWidth(gr->GetLineWidth()); 

  for (int i = 0; i < gr->GetN(); i++) 
  {
    toMollweide(gr->GetX()[i]*xscale, gr->GetY()[i], g.GetX()[i], g.GetY()[i],lon_0, reverseX); 
  }
}


void SkyMap::Paint(Option_t * opt) 
{

  if (range_set) 
  {
    double xmin, xmax, ymin,ymax; 

    toMollweide(range[1],range[3],xmax,ymax); 
    toMollweide(range[0],range[2],xmin,ymin); 

    if (xmax < xmin) std::swap(xmax,xmin); 

    sky_background.GetXaxis()->SetRangeUser(xmin,xmax);
    sky_background.GetYaxis()->SetRangeUser(ymin,ymax);
  }
  else
  {

    sky_background.GetXaxis()->UnZoom();
    sky_background.GetYaxis()->UnZoom(); 
    gPad->SetFrameLineColor(0); 
    gPad->SetFrameLineWidth(0); 
  }


  sky_background.Paint("COL A"); 

  paintSpecials(); 

  for (unsigned i = 0; i < sky_graphs.size(); i++) sky_graphs[i].Paint("lsame"); 
  for (unsigned i = 0; i < sky_markers.size(); i++) sky_markers[i].Paint("same"); 
  for (unsigned i = 0; i < sky_text.size(); i++) sky_text[i].Paint("same"); 
  left.Paint(); 
  right.Paint(); 

  if (zaxis && strcasestr(opt,"z") ) zaxis->Paint(); 

}


void SkyMap::setRange(double min_ra, double max_ra, double min_dec, double max_dec) 
{
  range_set = true; 
  range[0] = min_ra; 
  range[1] = max_ra; 
  range[2] = min_dec; 
  range[3] = max_dec; 
}


SkyMap::~SkyMap() 
{
  if (zaxis) delete zaxis; 
}


SkyMap::SkyMap(const SkyMap & other) 
  : TObject(other)
{
  reverseX = other.reverseX; 
  lon_0 = other.lon_0; 
  range_set = other.range_set; 
  if (range_set) 
  {
    std::copy(other.range, other.range+4, range); 
  }

  sky_background = other.sky_background; 
  sky_markers = other.sky_markers;
  sky_graphs = other.sky_graphs;
  for (const TLatex & ltx : other.sky_text)
  {
    sky_text.push_back(TLatex(ltx)); 
  }
  left = other.left; 
  right = other.right;
  if (other.zaxis) 
  {
    zaxis = new TPaletteAxis(*other.zaxis); 
  }
  else zaxis = 0; 
}
