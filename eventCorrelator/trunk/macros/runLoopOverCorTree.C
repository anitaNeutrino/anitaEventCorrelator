

void runLoopOverCorTree() {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("/usr/lib64/libfftw3.so");
  //  gSystem->Load("/unix/anita/softwareSLC4/install/lib/libfftw3.so");
  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaPlotter.so");
  gSystem->CompileMacro("loopOverCorTree.C","k");
  loopOverCorTree("/unix/anita1/rjn/peakCut5",1023);
  loopOverCorTree("/unix/anita1/rjn/peakCut5",1024);
  loopOverCorTree("/unix/anita1/rjn/peakCut5",1025);
  loopOverCorTree("/unix/anita1/rjn/peakCut5",1026);
  loopOverCorTree("/unix/anita1/rjn/peakCut5",1027);
  loopOverCorTree("/unix/anita1/rjn/peakCut5",1028);
}
