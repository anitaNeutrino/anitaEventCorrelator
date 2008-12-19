

void runPlotPrettyThings() {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/includes");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libfftw3.so");

  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaCorrelator.so");
  gSystem->Load("libRootFftwWrapper.so");
  gSystem->CompileMacro("plotPrettyThings.C","k");
  plotPrettyThingsEventNumber(3797,849,8); //Run,entry,antenna
  //  plotPrettyThings(1028,14289,16); //Run,entry,antenna
   //  plotPrettyThings(1028,33995,28);
}
