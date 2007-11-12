

void runCorTreeTest() {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  gSystem->AddIncludePath("-I${PLOTTER_DIR}");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  //  gSystem->Load("/usr/lib/libfftw3.so");
  gSystem->Load("/unix/anita/softwareSLC4/install/lib/libfftw3.so");
  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaPlotter.so");
  gSystem->CompileMacro("makeCorrelationTreeTest.C","k");
  makeCorrelationTreeTest(1028,146380,16); //Run,entry,antenna
   //  plotPrettyThings(1028,33995,28);
}
