

void runCorrTreeLoop1027() {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  gSystem->AddIncludePath("-I${PLOTTER_DIR}");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("/usr/lib64/libfftw3.so");
  //  gSystem->Load("/unix/anita/softwareSLC4/install/lib/libfftw3.so");
  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaCorrelator.so");
  gSystem->CompileMacro("correlationTreeLoop.C","k");
  correlationTreeLoop(1027);
  //   plotPrettyThings(1028,146380,16); //Run,entry,antenna
   //  plotPrettyThings(1028,33995,28);
}
