gSystem->Reset();

void runPlotEventLoc() {
  runPlotEventLocRemote(16,1267198,TMath::DegToRad()*35.1605,TMath::DegToRad()*226.987);
}

void runPlotEventLocRemote(int run, int event, double thetaWave, double phiWave) {
  //  gSystem->AddIncludePath(gSystem->ExpandPathName("-I${EVENT_READER_DIR}"));
  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");
  //  cout << gSystem->GetIncludePath() <<endl;
		
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libGeom.so");  
  gSystem->Load("libGraf3d.so");
  gSystem->Load("libRootFftwWrapper.so");     	  
  gSystem->Load("libAnitaEvent.so");   	  
  gSystem->Load("libAnitaCorrelator.so");   

  gSystem->CompileMacro("quickPlotEventLocation.C","k");
  quickPlotEventLocation("http://192.168.10.101/monitor/runs/",run,event,thetaWave,phiWave);

}
