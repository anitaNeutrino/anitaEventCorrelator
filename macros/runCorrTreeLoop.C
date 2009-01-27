

void runCorrTreeLoop() {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  //gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  //gSystem->AddIncludePath("-I${PLOTTER_DIR}");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  // gSystem->Load("/usr/lib64/libfftw3.so");
    gSystem->Load("/sw/lib/libfftw3.so");
     gSystem->Load("../../lib/libRootFftwWrapper.so");

    gSystem->Load("/Applications/root/lib/libMinuit.so");
  gSystem->Load("../../lib/libAnitaEvent.so");
  gSystem->Load("../../lib/libAnitaCorrelator.so");
  gSystem->CompileMacro("correlationTreeLoop.C","k");
   correlationTreeLoop(18);
  //  correlationTreeLoop(1028);
 

 //   plotPrettyThings(1028,146380,16); //Run,entry,antenna
   //  plotPrettyThings(1028,33995,28);
}
