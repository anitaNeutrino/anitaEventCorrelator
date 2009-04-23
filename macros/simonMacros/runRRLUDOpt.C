
void runRRLUDOpt() {

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

  gSystem->CompileMacro("rRLUDOpt.C","k");

  Double_t relDeltaOut=0;

  setAnt();
  setArray();

   TMinuit *myMin = new TMinuit(1);
   myMin->SetObjectFit(rRLUDOpt);
      myMin->SetFCN(iHateRoot);


       double startVal=0;               
       double stepSize=0.1;       
       double minVal=-2;           
       double maxVal=2;            
       Double_t p0 = 0;


       myMin->DefineParameter(0, "deltaR", startVal, stepSize, minVal, maxVal);
       myMin->DefineParameter(1, "deltaRL", startVal, stepSize, minVal, maxVal);
       myMin->DefineParameter(2, "deltaUD", startVal, stepSize, minVal, maxVal);
      
        myMin->FixParameter(1);  
       //myMin->FixParameter(2);  
 
       Double_t deltaR,deltaRErr;
       Double_t deltaRL,deltaRLErr;
       Double_t deltaUD,deltaUDErr;
  
//       //*********MINUIT METHOD*******************
       myMin->SetPrintLevel(-1);
       myMin->Migrad();   
       myMin->GetParameter(0,deltaR,deltaRLErr);
       myMin->GetParameter(1,deltaRL,deltaRLErr);
       myMin->GetParameter(2,deltaUD,deltaUDErr);

   
 	cout << deltaR <<"  "<< deltaRL << "  "<< deltaUD << endl;
      



  
}

void iHateRoot(Int_t& npar, Double_t* gin,
	       Double_t& f, Double_t* par, Int_t flag){


  double diffErr = rRLUDOpt(par);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}
