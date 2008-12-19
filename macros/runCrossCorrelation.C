
void runCrossCorrelation(int run,int entry) {

  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMinuit");
  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libRootFftwWrapper.so");
  gSystem->Load("libAnitaCorrelator.so");
  //  gSystem->Load("libRootFftwWrapper.so");
  gSystem->CompileMacro("/Users/Matt/WORK/eventCorrelator/trunk/macros/crossCorrelation.C","k");
  startCorrelation(run,entry); 
}
