

void runGroupDelay()
{
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");

    
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libGeom.so");  
  gSystem->Load("libMinuit.so");  
  gSystem->Load("libRootFftwWrapper.so");         
  gSystem->Load("libAnitaEvent.so");
  gSystem->CompileMacro("groupDelay.C","k");
  TChain *deltaChain = new TChain("deltaTTree");
  //  deltaChain->Add("/home/rjn/anita/data/deltaTTrees/justTaylorSvn220509/deltaTFile*.root");
  deltaChain->Add("/home/rjn/anita/data/deltaTTrees/justTaylorSimon290509WithDtHeading/deltaTFile*.root");
  groupDelay *t = new groupDelay(deltaChain);
  SetPaletteColours();
  t->Loop();



}
