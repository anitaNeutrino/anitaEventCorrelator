//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 



void runAngleOpt() {

  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");

  double startVal=0;               
  double stepSize=0.1;       
  double minVal=-2;           
  double maxVal=2;            
  Double_t p0 = 0;

  //Load libraries. Need to have ANITA_UTIL_INSTALL_DIR/lib and ROOTSYS/lib in the LD_LIBRARY_PATH               
  gSystem->Load("libfftw3.so");
  gSystem->Load("libMathMore.so");
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libGeom.so");  
  gSystem->Load("libMinuit.so");  
  gSystem->Load("libRootFftwWrapper.so");         
  gSystem->Load("libAnitaEvent.so");      
  gSystem->Load("libAnitaCorrelator.so");

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  gSystem->CompileMacro("anglePlotterOpt.C","k");
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(1);
  myMin->SetObjectFit(anglePlotterOpt);
  myMin->SetFCN(iHateRoot);


  for(int u = 0; u < 1; u++){

    int middle = 0;

      for(int y = 0; y <16; y++){

	int leftOpt, rightOpt;
	fGeomTool->getThetaPartners(middle,leftOpt,rightOpt); 


	myMin->DefineParameter(0, "antNum", middle, stepSize, minVal, maxVal);
	myMin->FixParameter(0);  
    
	myMin->DefineParameter(1, "deltaT", startVal, stepSize, minVal, maxVal);

  
	Double_t deltaT,deltaTErr;
  
	//*********MINUIT METHOD*******************
	myMin->SetPrintLevel(-1);
	myMin->Migrad();   
	myMin->GetParameter(1,deltaT,deltaTErr);

  

	setValue(rightOpt,deltaT);

	cout << middle << "  " << rightOpt << "  " << deltaT << endl;
  
	middle = rightOpt;

      }

  }

  
  //   myMin->DeleteArrays();

  //   myMin->DeleteArrays();

  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){


  double diffErr = anglePlotterOpt(par);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


// 0  9  -0.00415702
// 9  1  0.0442874
// 1  10  -0.249513
// 10  2  0.493758
// 2  11  -0.550521
// 11  3  0.501072
// 3  12  -0.732701
// 12  4  0.993376
// 4  13  -0.95656
// 13  5  0.862654
// 5  14  -0.911761
// 14  6  1.02809
// 6  15  -0.982902
// 15  7  0.932202
// 7  8  -1.12693



