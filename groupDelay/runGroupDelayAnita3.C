

void runGroupDelayAnita3()
{
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");

    
  //  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libGeom.so");  
  gSystem->Load("libMinuit.so");  
  gSystem->Load("libRootFftwWrapper.so");         
  gSystem->Load("libAnitaEvent.so");
  gSystem->CompileMacro("groupDelayAnita3.C","k");

  TFile *finput = new TFile("offsetdelay_WIDE_VPOL.root", "read");
  //  TFile *finput = new TFile("offsetdelay_WIDE.root", "read");
  TH2D *histDtPhi = (TH2D*)finput->Get("deltaTvsPhiDiff");
  histDtPhi->SetName("histDtPhi");

  char histName[100];
  TH2D *histDtPhiAnt[48];
   for(int ant=0;ant<48;ant++) {
     sprintf(histName,"histDtPhiAnt%d",ant);
     histDtPhiAnt[ant]= (TH2D*)finput->Get(histName);
   }

   groupDelayAnita3(histDtPhi, histDtPhiAnt);




}
