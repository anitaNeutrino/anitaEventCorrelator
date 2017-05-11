#include "AntarcticaGeometry.h" 
#include "AnitaGeomTool.h" 
#include "RampdemReader.h" 

#include "TRandom.h" 
#include "TGraph2D.h" 
#include "TString.h" 

void AntarcticCoord::asString(TString * s) const
{

  if (type == CARTESIAN) 
  {
    s->Form("%g N %g E %g m",x,y,z); 
  }
  
  if (type == STEREOGRAPHIC)
  {
    s->Form("%g m(E), %g m(N), %g m",x,y,z); 
  }

  if (type == CARTESIAN)
  {
    s->Form("(%g,%g,%g) m",x,y,z); 
  }

}

void AntarcticCoord::convert(CoordType t)
{

  if (type == WGS84) 
  {
    if (t == CARTESIAN) 
    {
      double cartesian[3]; 
      AnitaGeomTool::Instance()->getCartesianCoords(x,y,z, cartesian); 
      x = cartesian[0]; 
      y = cartesian[1]; 
      z = cartesian[2]; 
    }

    else if ( t==STEREOGRAPHIC) 
    {
      RampdemReader::LonLatToEastingNorthing(y,x,x,y); 
    }

  }

  else if (type == STEREOGRAPHIC)
  {


    if (t == WGS84) 
    {
      RampdemReader::EastingNorthingToLonLat(x,y,y,x); 
    }

    else if (t == CARTESIAN) 
    {
      //be lazy 
      convert(WGS84); 
      convert(CARTESIAN); 
    }

  }

  else if (type == CARTESIAN) 
  {

    if (t == WGS84) 
    {
      double cartesian[3]; 
      cartesian[0] = x; 
      cartesian[1] = y; 
      cartesian[2] = z; 
      AnitaGeomTool::Instance()->getLatLonAltFromCartesian(cartesian, x,y,z); 
    }
    else if ( t== STEREOGRAPHIC)
    {
      //be lazy 
      convert(WGS84); 
      convert(STEREOGRAPHIC); 
    }
  }

  else 
  {
    fprintf(stderr,"Invalid CoordType: %d!", type); 
    return; 
  }


  type = t; 

}



AntarcticSegmentationScheme * AntarcticSegmentationScheme::factory(const char * key) 
{

  if (strstr(key,"stereographic_grid"))
  {
    int nx,ny,maxE,maxN; 
    ny = -1; 
    nx = 0; 
    maxE = 3330000; 
    maxN = 3330000; 
    sscanf(key,"stereographic_grid_%d_%d_%d_%d", &nx, &ny, &maxN, &maxE); 
    if (!nx) return 0; 
    if (ny < 0) ny = nx; 
    return new StereographicGrid(nx,ny,maxE,maxN); 
  }

  return 0; 
}

const int nsamples = 64; 

void AntarcticSegmentationScheme::Draw(const char * opt) const 
{
  int N = NSegments(); 

  TGraph2D * g  = new TGraph2D(N*nsamples); 
  g->SetBit(kCanDelete); 

  AntarcticCoord samples[nsamples]; 

  for (int i =0; i < N; i++) 
  {
    /* do not need altitude for this! */ 
    sampleSegment(i, nsamples, samples, true,false); 

    for (int j =0;  j<  nsamples; j++) 
    {
      samples[j].to(AntarcticCoord::STEREOGRAPHIC); 
      g->SetPoint(nsamples*i+j,samples[j].x,samples[j].y,i); 
    }
  }

  g->Draw(opt); 
}



StereographicGrid::StereographicGrid(int NX, int NY, double MAX_E, double MAX_N)
  : nx(NX), ny(NY), max_E(MAX_E), max_N(MAX_N) 
{

  dx = (2*max_E)/nx; 
  dy = (2*max_N)/ny; 
}

int StereographicGrid::getSegmentIndex(const AntarcticCoord & coord) const
{
  AntarcticCoord stereo = coord.as(AntarcticCoord::STEREOGRAPHIC); 

  if (stereo.x >  max_E || stereo.x < -max_E || stereo.y > max_N || stereo.y < -max_N) return -1; 


  int xbin = nx/2 + stereo.x / dx; 
  int ybin = ny/2 - stereo.y / dy; 

  return xbin + nx * ybin; 

}

void StereographicGrid::getSegmentCenter(int idx, AntarcticCoord * fillme, bool fillalt) const 
{

  int xbin = idx % nx; 
  int ybin = idx / nx; 
  double x = (xbin +0.5) *dx - max_E; 
  double y = max_N - (ybin +0.5) *dy ; 
  double z = fillalt ? RampdemReader::SurfaceAboveGeoidEN(x,y,dataset): 0; 
  fillme->set(AntarcticCoord::STEREOGRAPHIC,x,y,z); 
}


AntarcticCoord * StereographicGrid::sampleSegment(int idx, int N, AntarcticCoord * fill, bool random, bool fillalt) const 
{

  if (!fill) fill = new AntarcticCoord[N]; 

  int xbin = idx % nx; 
  int ybin = idx / nx; 
  double lx = (xbin) *dx - max_E; 
  double ly = max_N - (ybin) *dy ; 


  if (random) 
  {
    for (int i = 0; i < N; i++)
    {
      double x = gRandom->Uniform(lx, lx+dx); 
      double y = gRandom->Uniform(ly, ly+dy); 
      double z = fillalt ? RampdemReader::SurfaceAboveGeoidEN(x,y,dataset): 0; 
      fill[i].set(AntarcticCoord::STEREOGRAPHIC,x,y,z); 
    }
  }
  else
  {
    int grid = sqrt(N) + 0.5 ; 

    for (int i = 0; i < N; i++)
    {
      //TODO check if this is what i want! 
      double x = lx + (dx / (grid+1)) * (1+(i % grid)); 
      double y = ly + (dy / (grid+1)) * (1+(i / grid)); 
      double z = fillalt ? RampdemReader::SurfaceAboveGeoidEN(x,y,dataset): 0; 
      fill[i].set(AntarcticCoord::STEREOGRAPHIC,x,y,z); 
    }
  }

  return fill; 
}



void StereographicGrid::Draw(const char * opt) const
{
  TH2I h("tmp","Stereographic Grid", nx, -max_E, max_E, ny, -max_N, max_N); 
  for (int i = 1; i <= nx; i++) 
  {
    for (int j = 1; j <= ny; j++) 
    {
      AntarcticCoord c(AntarcticCoord::STEREOGRAPHIC, h.GetXaxis()->GetBinCenter(i), h.GetYaxis()->GetBinCenter(j));
      h.SetBinContent(i,j,  getSegmentIndex(c)); 
    }
  }
  h.DrawCopy(opt); 
}



