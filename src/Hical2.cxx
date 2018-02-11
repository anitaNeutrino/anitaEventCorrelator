#include "Hical2.h"

using namespace std;

ClassImp(Hical2)


Hical2::Hical2(){}

double Hical2::isHical(UInt_t eventNumber, int geomCut){
  //initialize the hical stuff
  if(INIT_HICAL==0){
    initHical();
    INIT_HICAL=1;
  }
  //quick sanity checks for hical not being in the air
  if(eventNumber<31059701){
    //hical is not not in the air yet!
    return 0;
  }
  if(eventNumber>72073425){
    // hical is dead!;
    return 0;
  }
  //quick check to see if both payloads are off
  adset->getEvent(eventNumber);
  if(hc2aOn(adset->header()->triggerTime)==0&&hc2bOn(adset->header()->triggerTime)==0){
    //    cout<<"both are off"<<endl;
    return 0;
  }

  if(geomCut==1){
    //variables for the pointing to HC
    // double a4_hc2a, a4_hc2b;  
    // whereAreHical(eventNumber, &a4_hc2a, &a4_hc2b);
    // cout<<a4_hc2a<<" "<<a4_hc2b<<" "<<FFTtools::wrap(sum->anitaLocation.heading-sum->peak[0][0].phi, 360)<<endl;
    
    // if(abs(a4_hc2a-FFTtools::wrap(sum->anitaLocation.heading-sum->peak[0][0].phi, 360))<2.5){
    
    //   return 1;
    // }
    // else if(abs(a4_hc2b-FFTtools::wrap(sum->anitaLocation.heading-sum->peak[0][0].phi, 360))<2.5){
    
    //   return 1;
    // }
    // else {
    //   return 0;
    // }
    return 0;
  }

  return -1;
}
double Hical2::isHical(AnitaEventSummary *sum, int geomCut){
  //quick sanity checks for hical not being in the air
  if(sum->eventNumber<31059701){
    //hical is not not in the air yet!
    return 0;
  }
  if(sum->eventNumber>72073425){
    // hical is dead!;
    return 0;
  }
  //quick check to see if both payloads are off
  //  ad->getEvent(sum->eventNumber);
  if(hc2aOn(sum->realTime)==0&&hc2bOn(sum->realTime)==0){
    //    cout<<"both are off"<<endl;
    return 0;
  }

  if(geomCut==1){
    //variables for the pointing to HC
    double a4_hc2a, a4_hc2b;  
    whereAreHical(sum->eventNumber, &a4_hc2a, &a4_hc2b);
    //cout<<a4_hc2a<<" "<<a4_hc2b<<" "<<FFTtools::wrap(sum->anitaLocation.heading-sum->peak[0][0].phi, 360)<<endl;
    
    if(abs(a4_hc2a-FFTtools::wrap(sum->anitaLocation.heading-sum->peak[0][0].phi, 360))<2.5&&hc2aOn(sum->realTime)==true){
    
      return 1;
    }
    else if(abs(a4_hc2b-FFTtools::wrap(sum->anitaLocation.heading-sum->peak[0][0].phi, 360))<2.5&&hc2bOn(sum->realTime)==true){
    
      return 1;
    }
    else {
      return 0;
    }
  }

  return -1;
}

int Hical2::whereAreHical(UInt_t eventNumber, double * angleToA, double * angleToB){
  if(eventNumber<31059701){
    //hical is not not in the air yet!
    *angleToA=-999;
    *angleToB=-999;
    return 0;
  }
  if(eventNumber>72073425){
    *angleToA=-999;
    *angleToB=-999;
    // hical is dead!;
    return 0;
  }

  static AnitaDataset *ad=0;
  if(!ad)ad=new AnitaDataset(157);
  
  static TFile *hc2ahk_file = 0;
  static TTree *hc2ahk_tree = 0;
  static HCHKTree *hc2ahk=new HCHKTree();
  if(!hc2ahk_file){
    char filename[1000];
    char *dir=getenv("ANITA_UTIL_INSTALL_DIR");
    sprintf(filename,"%s/share/anitaCalib/hc2ahk.root",dir);
    hc2ahk_file=TFile::Open(filename);
    hc2ahk_tree = (TTree*)hc2ahk_file->Get("tree");
    hc2ahk_tree->SetBranchAddress("hktree", &hc2ahk);
    hc2ahk_tree->BuildIndex("hktree.time");
  }

  static TFile *hc2bhk_file = 0;
  static TTree *hc2bhk_tree = 0;
  static HCHKTree *hc2bhk = new HCHKTree();
  if(!hc2bhk_file){
    char filename[1000];
    char *dir=getenv("ANITA_UTIL_INSTALL_DIR");
    sprintf(filename,"%s/share/anitaCalib/hc2bhk.root",dir);
    hc2bhk_file=TFile::Open(filename);
    hc2bhk_tree = (TTree*)hc2bhk_file->Get("tree");
    hc2bhk_tree->SetBranchAddress("hktree", &hc2bhk);
    hc2bhk_tree->BuildIndex("hktree.time");
  }

 
  
  ad->getEvent(eventNumber);
  double atime1, alat1, alon1, aalt1, atime2, alat2, alon2, aalt2;
  double btime1, blat1, blon1, balt1, btime2, blat2, blon2, balt2;
  double alat, alon, aalt, blat, blon, balt;
  double frachk;
  
  //for hc2a
  int aentry = hc2ahk_tree->GetEntryNumberWithBestIndex(ad->header()->triggerTime, ad->header()->triggerTime);
  hc2ahk_tree->GetEntry(aentry);
  if(aentry==hc2ahk_tree->GetEntries()-1){
    *angleToA=-999;
  }
  else{
    atime1 = hc2ahk->time;
    alat1  = hc2ahk->lat;
    alon1  = hc2ahk->lon;
    aalt1  = hc2ahk->alt;

    hc2ahk_tree->GetEntry(aentry+1);

    atime2 = hc2ahk->time;
    alat2  = hc2ahk->lat;
    alon2  = hc2ahk->lon;
    aalt2  = hc2ahk->alt;

    frachk = ((double)ad->header()->triggerTime+((double)ad->header()->triggerTimeNs/1000000000.)-atime1)/(atime2-atime1);

    alat = (frachk*(alat2-alat1))+alat1;
    alon = (frachk*(alon2-alon1))+alon1;
    aalt = (frachk*(aalt2-aalt1))+aalt1;
  
    *angleToA = angleToThing(ad->gps()->latitude, ad->gps()->longitude, alat, alon);
  }
  //for hc2b
  int bentry = hc2bhk_tree->GetEntryNumberWithBestIndex(ad->header()->triggerTime, ad->header()->triggerTime);
  hc2bhk_tree->GetEntry(bentry);

  if(bentry==hc2bhk_tree->GetEntries()-1){
    *angleToB=-999;
  }
  else{
    btime1 = hc2bhk->time;
    blat1  = hc2bhk->lat;
    blon1  = hc2bhk->lon;
    balt1  = hc2bhk->alt;

    hc2bhk_tree->GetEntry(bentry+1);

    btime2 = hc2bhk->time;
    blat2  = hc2bhk->lat;
    blon2  = hc2bhk->lon;
    balt2  = hc2bhk->alt;

    frachk = ((double)ad->header()->triggerTime+((double)ad->header()->triggerTimeNs/1000000000.)-btime1)/(btime2-btime1);

    blat = (frachk*(blat2-blat1))+blat1;
    blon = (frachk*(blon2-blon1))+blon1;
    balt = (frachk*(balt2-balt1))+balt1;
  
    *angleToB = angleToThing(ad->gps()->latitude, ad->gps()->longitude, blat, blon);
  }
  
  return 1;
}

double Hical2::dPhi(int aorb, int peak){
  return .5;
}

int Hical2::initHical(){
  AnitaVersion::set(4);
  adset = new AnitaDataset(157);
  static TTree *hc2ahk_tree = 0;
  static HCHKTree *hc2ahk=new HCHKTree();

  static TTree *hc2bhk_tree = 0;
  static HCHKTree *hc2bhk=new HCHKTree();

  char filename[1000];
  char *dir=getenv("ANITA_UTIL_INSTALL_DIR");

  sprintf(filename,"%s/share/anitaCalib/hc2ahk.root",dir);
  TFile *hc2ahk_file = TFile::Open(filename);
  hc2ahk_tree = (TTree*)hc2ahk_file->Get("tree");
  hc2ahk_tree->SetBranchAddress("hktree", &hc2ahk);
  hc2ahk_tree->BuildIndex("hktree.time");

  sprintf(filename, "%s/share/anitaCalib/hc2bhk.root",dir);
  TFile *hc2bhk_file = TFile::Open(filename);
  hc2bhk_tree = (TTree*)hc2bhk_file->Get("tree");
  hc2bhk_tree->SetBranchAddress("hktree", &hc2bhk);
  hc2bhk_tree->BuildIndex("hktree.time");
  
  return 1;
}


double Hical2::angleToThing(double a4lat, double a4lon, double hclat, double hclon){
  double lat1=deg2rad(a4lat);
  double lat2=deg2rad(hclat);
  double lon1 = deg2rad(a4lon);
  double lon2 = deg2rad(hclon);
  double x = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon2-lon1);
  double y = sin(lon2-lon1) * cos(lat2);
  double bearing = atan2(y, x);
  if(bearing<0)bearing=((2.*TMath::Pi())+bearing);
  return rad2deg(bearing);
}

double Hical2::deg2rad(double deg) {
  return (deg * TMath::Pi() / 180.);
}


double Hical2::rad2deg(double rad) {
  return (rad * 180. / TMath::Pi());
}

//hc2a runs, known times that the payload is on or off from multiple sources
bool Hical2::hc2aOn(UInt_t triggerTime){
  static TFile *infile = 0;
  if(!infile){
    char filename[1000];
    char *dir=getenv("ANITA_UTIL_INSTALL_DIR");
    sprintf(filename,"%s/share/anitaCalib/hc2aruns.root",dir);
    infile = TFile::Open(filename);
  }
  //static TTree * intree = 0;
  UInt_t ontime, offtime;
  bool on=false;
  // if(!intree){
   TTree* intree=(TTree*)infile->Get("tree");
    
    intree->SetBranchAddress("ontime", &ontime);
    intree->SetBranchAddress("offtime", &offtime);
    //  }
  for(int i=0;i<intree->GetEntries();i++){
    intree->GetEntry(i);
    if(triggerTime>=ontime&&triggerTime<=offtime){
      on=true;
    }
  }

  UInt_t t = triggerTime-1480000000;
  if((t>717055&&t<718593)||
     (t>787583&&t<791809)||
     (t>1572095&&t<1627009)||
     (t>1705855&&t<1705985)||
     (t>1748991&&t<1750529)||
     (t>1775615&&t<1776897)||
     (t>1811199&&t<1813633)||
     (t>1861375&&t<1868289)||
     (t>1912319&&t<1914241)||
     (t>1951615&&t<1953537)||
     (t>1988991&&t<1989121)||
     (t>2012287&&t<2015617)||
     (t>2042111&&t<2043137))
    on=true;  
  //infile->Close();
  return on;
}

//hical 2b
bool Hical2::hc2bOn(UInt_t triggerTime){
  static TFile *infile = 0;
  if(!infile){
    char filename[1000];
    char *dir=getenv("ANITA_UTIL_INSTALL_DIR");
    sprintf(filename,"%s/share/anitaCalib/hc2bruns.root",dir);
    infile = TFile::Open(filename);
  }
  //  static TTree * intree = 0;
  UInt_t ontime, offtime;
  bool on=false;
  //if(!intree){
    TTree * intree=(TTree*)infile->Get("tree");
    intree->SetBranchAddress("ontime", &ontime);
    intree->SetBranchAddress("offtime", &offtime);
    //  }
  for(int i=0;i<intree->GetEntries();i++){
    intree->GetEntry(i);
    if(triggerTime>=ontime&&triggerTime<=offtime){
      on=true;
    }
  }

    UInt_t t = triggerTime-1480000000;
  if((t>1494143&&t<1514753)||
     (t>1564159&&t<1612801)||
     (t>1646079&&t<1648129)||
     (t>1704703&&t<1705985)||
     (t>1724159&&t<1725697)||
     (t>1749503&&t<1750273)||
     (t>1776255&&t<1777281)||
     (t>1811199&&t<1813633)||
     (t>1861375&&t<1867905)||
     (t>1912319&&t<1914113)||
     (t>1951743&&t<1953409)||
     (t>1989631&&t<1990785)||
     (t>2041855&&t<2043265)||
     (t>2053375&&t<2054657)||
     (t>2134911&&t<2138113)||
     (t>2162175&&t<2164225)||
     (t>2177919&&t<2181249)||
     (t>2214399&&t<2216321)||
     (t>2251391&&t<2253185)||
     (t>2265087&&t<2267265))
    on=true;

  //infile->Close();
  return on;
}
