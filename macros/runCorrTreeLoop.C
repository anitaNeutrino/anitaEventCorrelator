

void runCorrTreeLoop() {
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");

  //Load libraries. Need to have ANITA_UTIL_INSTALL_DIR/lib and ROOTSYS/lib in the LD_LIBRARY_PATH               
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libGeom.so");  
  gSystem->Load("libMinuit.so");  
  gSystem->Load("libRootFftwWrapper.so");         
  gSystem->Load("libAnitaEvent.so");      
  gSystem->Load("libAnitaCorrelator.so");


  gSystem->CompileMacro("correlationTreeLoop.C","k");
  correlationTreeLoop(18,"http://www.hep.ucl.ac.uk/uhen/anita/private/monitor2/runs/fromLoki/","/home/rjn/anita/data/corTrees/","/home/rjn/anita/data/deltaTTrees");
}
