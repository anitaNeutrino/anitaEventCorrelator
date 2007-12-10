

void runLoopOverInitialCorTree() {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${EVENT_READER_DIR}");
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/includes");
  //  cout << gSystem->GetIncludePath() <<endl;
			  
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("/usr/lib64/libfftw3.so");
  //  gSystem->Load("/unix/anita/softwareSLC4/install/lib/libfftw3.so");
  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaPlotter.so");
  gSystem->CompileMacro("loopOverCorTree.C","k");
  char fileName[180];
  FileStat_t staty;
  for(int run=1030;run<=1239;run++) {
     sprintf(fileName,"/unix/anita1/rjn/initialPeakCut5/corRun%d.root",run);
     if(gSystem->GetPathInfo(fileName,staty)) {
	continue;
	//cout << fileName << endl;
     }
     loopOverCorTree("/unix/anita1/rjn/initialPeakCut5",run);
  }
}
