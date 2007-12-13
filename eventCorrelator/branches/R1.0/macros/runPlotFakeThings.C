

void runPlotFakeThings() {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  gSystem->AddIncludePath("-I${PLOTTER_DIR}");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("/usr/lib64/libfftw3.so");
  //  gSystem->Load("/unix/anita/softwareSLC4/install/lib/libfftw3.so");
  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaPlotter.so");
  gSystem->CompileMacro("plotFakePrettyThings.C","k");
  //  plotFakePrettyThings(0,16); //entry,antenna  
  plotFakePrettyThings(0,23); //entry,antenna

}
