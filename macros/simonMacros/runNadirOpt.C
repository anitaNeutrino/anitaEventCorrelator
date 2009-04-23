//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 



void runNadirOpt() {

  gSystem->AddIncludePath("-I${ANITA_UTIL_INSTALL_DIR}/include");

  double startVal=0;               
  double stepSize=0.1;       
  double minVal=-0.5;           
  double maxVal=0.5;            
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
  gSystem->CompileMacro("nadirOpt.C","k");
  
  Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(150);
  myMin->SetObjectFit(nadirOpt);
  myMin->SetFCN(iHateRoot);
  //setArray();



 for(int y = 0; y <8; y++){

   char name[30];
   sprintf(name,"r%d",y);
	myMin->DefineParameter(y, name, startVal, stepSize, minVal, maxVal);
   sprintf(name,"z%d",y);
	myMin->DefineParameter(y+8, name, startVal, stepSize, minVal, maxVal);
	   sprintf(name,"phi%d",y);
	myMin->DefineParameter(y+16, name, startVal, stepSize, minVal, maxVal);
 }


  
	Double_t deltaR[8],deltaRErr[8];
	Double_t deltaZ[8],deltaZErr[8];
	Double_t deltaPhi[8],deltaPhiErr[8];
  
	//*********MINUIT METHOD*******************
	myMin->SetPrintLevel(-1);
	myMin->Migrad();   

// 	for(int u = 0; u <8; u++){
// 	  myMin->GetParameter(u,deltaR[u],deltaRErr[u]);
// 	cout << "deltaR[" << u << "] = " << deltaR[u] << ";" << endl;
// 	}
// 	for(int u = 0; u <8; u++){
// 	  myMin->GetParameter(u+8,deltaZ[u],deltaZErr[u]);
// 	cout << "deltaZ[" << u << "] = " << deltaZ[u] << ";" << endl;
// 	}
// 	for(int u = 0; u <8; u++){
// 	  myMin->GetParameter(u+16,deltaPhi[u],deltaPhiErr[u]);
// 	cout << "deltaPhi[" << u << "] = " << deltaPhi[u] << ";" << endl;
// 	}


	for(int u = 0; u <8; u++){
	  myMin->GetParameter(u,deltaR[u],deltaRErr[u]);
	  //cout << "deltaR[" << u << "] = " << deltaR[u] ;
	
	
	  myMin->GetParameter(u+8,deltaZ[u],deltaZErr[u]);
	  //cout << " deltaZ[" << u << "] = " << deltaZ[u] ;
	
	  myMin->GetParameter(u+32,deltaPhi[u],deltaPhiErr[u]);
	  //cout << " deltaPhi[" << u << "] = " << deltaPhi[u] << ";" << endl;

	  cout << u << "  " << deltaPhi[u]<< "  " << deltaR[u]<< "  " << deltaZ[u]<< "  " << 0 << endl;


	}

  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){


  double diffErr = nadirOpt(par);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


