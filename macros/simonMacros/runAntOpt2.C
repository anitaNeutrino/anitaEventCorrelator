


void runAntOpt2() {

  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");

  double startVal=0;               
  double stepSize=0.1;       
  double minVal=-2;           
  double maxVal=2;            
  Double_t p0 = 0;

  Double_t deltaR[16],deltaRErr[16];
  // Double_t relDelta[16],relDeltaErr[16];
  
  Double_t deltaZ[16],deltaZErr[16];
  //Double_t relDeltaZ[16],relDeltaZErr[16];
  
  Double_t deltaPhi[16],deltaPhiErr[16];

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
  gSystem->CompileMacro("antOpt2.C","k");

  Double_t relDeltaOut=0;

  setAnt();
  setArray();

  TMinuit *myMin = new TMinuit(1);
  myMin->SetObjectFit(antOpt2);
  myMin->SetFCN(iHateRoot);


  for(int u = 0; u < 1; u++){

    int middle = 0;

       for(int y = 0; y <1; y++){


     

      int leftOpt, rightOpt;
      fGeomTool->getThetaPartners(middle,leftOpt,rightOpt); 

      double startVal=0;               
      double stepSize=0.1;       
      double minVal=-2;           
      double maxVal=2;            
      Double_t p0 = 0;

      char deltaRPar[100];
      char deltaZPar[100];
      char deltaPhiPar[100];

      myMin->DefineParameter(0, "antNum", middle, stepSize, minVal, maxVal);
      myMin->FixParameter(0);  

//       for(int r = 0; r < 16; r++){
// 	sprintf(deltaRPar,"deltaR%d",r);
// 	myMin->DefineParameter(r+5, deltaRPar, startVal, stepSize, minVal, maxVal);
//       }

//       for(int z = 0; z < 16; z++){
// 	sprintf(deltaZPar,"deltaZ%d",z);
// 	// myMin->DefineParameter(z+16, deltaZPar, startVal, stepSize, minVal, maxVal);
// 	//myMin->FixParameter(z+16);  
//       }
  
      for(int phi = 0; phi < 16; phi++){
	sprintf(deltaPhiPar,"deltaPhi%d",phi);
	myMin->DefineParameter(phi+5, deltaPhiPar, startVal, stepSize, minVal, maxVal);
	//myMin->FixParameter(phi+1);  
      }

  
    
      myMin->DefineParameter(1, "deltaR", startVal, stepSize, minVal, maxVal);
      myMin->DefineParameter(2, "deltaZ", startVal, stepSize, minVal, maxVal);
      myMin->DefineParameter(3, "deltaPhi", startVal, stepSize, minVal, maxVal);
      myMin->DefineParameter(4, "deltaT", startVal, stepSize, minVal, maxVal);
      myMin->FixParameter(4);
      myMin->FixParameter(1);
      myMin->FixParameter(2);
      myMin->FixParameter(3);

 
      Double_t deltaT,deltaTErr;
      Double_t deltaRO,deltaRErrO;
      Double_t deltaZO,deltaZErrO;
      Double_t deltaPhiO,deltaPhiErrO;
  
      //*********MINUIT METHOD*******************
      myMin->SetPrintLevel(-1);
      myMin->Migrad();   
      myMin->GetParameter(1,deltaRO,deltaRErrO);
      myMin->GetParameter(2,deltaZO,deltaZErrO);
      myMin->GetParameter(3,deltaPhiO,deltaPhiErrO);
      myMin->GetParameter(4,deltaT,deltaTErr);
  
      for(int r = 0; r < 16; r++){
	myMin->GetParameter(r+5,deltaR[r],deltaRErr[r]);
      }
   
//       for(int z = 0; z < 16; z++){
// 	// myMin->GetParameter(z+16,deltaZ[z],deltaZErr[z]);
//       }

//       for(int phi = 0; phi < 16; phi++){
// 	myMin->GetParameter(phi+16,deltaPhi[phi],deltaRErr[phi]);
//       }


      for(int r = 0; r < 16; r++){
	cout << deltaR[r] <<"  "<< 0 << "  "<< 0 << endl;
      }  


      setAntValue(rightOpt,deltaRO,deltaPhiO,deltaZO);
	setValue(rightOpt,deltaT);



      //	cout << middle << "  " << rightOpt << "  " << deltaT << endl;
  
      //	cout << "deltaTArrayMod[" << rightOpt << "] = " << deltaT << ";" << endl;

      middle = rightOpt;

    }

    //  	printAntArray();
    //	cout << "  " << endl;
    // 	printArray();

  }

  
  //   myMin->DeleteArrays();

  //   myMin->DeleteArrays();

  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){


  double diffErr = antOpt2(par);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}
