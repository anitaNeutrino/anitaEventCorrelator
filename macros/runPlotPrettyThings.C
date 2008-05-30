

void runPlotPrettyThings() {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("/usr/lib64/libfftw3.so");

  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaPlotter.so");
  gSystem->Load("libRootFftwWrapper.so");
  gSystem->CompileMacro("plotPrettyThings.C","k");
  plotPrettyThingsEventNumber(1032,612849,8); //Run,entry,antenna
  //  plotPrettyThings(1028,14289,16); //Run,entry,antenna
   //  plotPrettyThings(1028,33995,28);
}
