
//TTree *headTree; 
//TTree *adu5PatTree; 
//TTree *corTree; 

void runAntOpt(Int_t antIn1, Int_t antIn2) {
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


   gSystem->CompileMacro("antOpt.C","k");
   gSystem->CompileMacro("correlationTreeLoopOpt.C","k");
  
  ifstream optFile("optFile.txt");
     
  if (optFile.fail()) {
    cout << "optFile.txt" << "  file  NOT FOUND" << endl;
    return 1;
  }
  
  else{
    
  Double_t deltaROut=0;    
  Double_t deltaZOut=0;    
  Double_t deltaPhiOut=0;    
  optFile >> deltaROut;
  optFile >> deltaZOut;
  optFile >> deltaPhiOut;

  optFile.close();
  ofstream optFileOut("optFile.txt");

   Double_t relDeltaOut=0;

  TMinuit *myMin = new TMinuit(1);
  myMin->SetObjectFit(antOpt);
  myMin->SetFCN(iHateRoot);



  //myMin->DefineParameter(0, "deltaR", startVal, stepSize, minVal, maxVal);
  //myMin->DefineParameter(1, "relDelta", deltaROut, stepSize, minVal, maxVal);
  

  //myMin->DefineParameter(3, "deltaZ", startVal, stepSize, minVal, maxVal);
  //myMin->DefineParameter(4, "relDeltaZ", deltaZOut, stepSize, minVal, maxVal);
  
  //myMin->DefineParameter(5, "deltaPhi", startVal, stepSize, minVal, maxVal);
  //myMin->DefineParameter(6, "relDeltaPhi", deltaPhiOut, stepSize, minVal, maxVal);

 char deltaRPar[100];
 char deltaZPar[100];
 char deltaPhiPar[100];

  for(int r = 0; r < 1; r++){
    sprintf(deltaRPar,"deltaR%d",r);
    myMin->DefineParameter(r, deltaRPar, startVal, stepSize, minVal, maxVal);
  }

  for(int z = 0; z < 16; z++){
    sprintf(deltaZPar,"deltaZ%d",z);
    // myMin->DefineParameter(z+16, deltaZPar, startVal, stepSize, minVal, maxVal);
    //myMin->FixParameter(z+16);  
  }
  
  for(int phi = 0; phi < 16; phi++){
    sprintf(deltaPhiPar,"deltaPhi%d",phi);
    myMin->DefineParameter(phi+1, deltaPhiPar, startVal, stepSize, minVal, maxVal);
    myMin->FixParameter(phi+1);  
  }
 
  myMin->DefineParameter(18, "deltaHeading", startVal, stepSize, minVal, maxVal);

  //myMin->DefineParameter(33, "ant", antIn1, stepSize, antIn1-1, antIn1+1);

 // myMin->FixParameter(1);  
 // myMin->FixParameter(2);  
 // myMin->FixParameter(4);  
 // myMin->FixParameter(6);  

  Double_t deltaR[16],deltaRErr[16];
  // Double_t relDelta[16],relDeltaErr[16];
  
  Double_t deltaZ[16],deltaZErr[16];
  //Double_t relDeltaZ[16],relDeltaZErr[16];
  
   Double_t deltaPhi[16],deltaPhiErr[16];
   //Double_t relDeltaPhi[16],relDeltaPhiErr[16];
   
Double_t deltaHeading,deltaHeading;

   //*************VARY R METHOD*******************

//    Double_t par2[32];

//    //TH1F *rGraph = new TH1F("rGraph", "rGraph", 201,-0.2,0.2);
//    int count = 0;
//    vector<float> diffErrVec;
//    vector<float> rDeltaVec;
//      for(float rDelta = -0.5; rDelta<0.5; rDelta = rDelta+0.1){
//        //cout << rDelta << endl;
//        for(int r = 0; r < 1; r++){ 
// 	 par2[r] = rDelta;
//        }
       
//        for(int phi = 0; phi < 16; phi++){
// 	 // par2[phi+16]=0;
//        }
       
//        float diffErr = (float)antOpt(par2));
//      diffErrVec.push_back(diffErr);
//      rDeltaVec.push_back(rDelta);
//        count++;
//        //rGraph->Fill(rDelta,diffErr);
//      }

//   Float_t diffErrArray[40];
//   Float_t rDeltaArray[40];
//   for(int i = 0; i < count; i++){
//      diffErrArray[i] = diffErrVec[i];
//      rDeltaArray[i] = rDeltaVec[i];
//   }

//   TCanvas *canHead = new TCanvas("canHead","canHead",600,400);
//   TGraph *rGraph = new TGraph(count,rDeltaArray,diffErrArray);
//   rGraph->Draw("ap");
//   rGraph->GetXaxis()->SetTitle("delta R (m)");
//   rGraph->GetYaxis()->SetTitle("Sum of (actual - expected) time (s)");
//   rGraph->SetTitle("Calculation of the sum of actual - expected time for all pairs");
//   rGraph->SetMarkerStyle(22);
//   rGraph->SetMarkerColor(1);
//   TF1 *f1 = new TF1("f1","x*[0]+[1]",-0.5,0.5);
//   f1->SetLineWidth(0.1);  
//   rGraph->Fit("f1","R");
//   canHead->Update();

//*********VARY R PHASE CENTRE *****************

//  Float_t diffErrArray[9];
//  Float_t rDeltaArray[9];
//   Float_t rmsDeltaArray[9];

//   diffErrArray[0] = -0.2;    
//   rDeltaArray[0] = -0.54523;    
//   rmsDeltaArray[0] = 1.61264;

//   diffErrArray[1] = -0.15;  
//   rDeltaArray[1] = -0.449344;   
//   rmsDeltaArray[1] = 1.42441;

//   diffErrArray[2] = -0.1;   
//   rDeltaArray[2] =  -0.352232;  
//   rmsDeltaArray[2] =  1.32794;

//   diffErrArray[3] =  -0.05;  
//   rDeltaArray[3] = -0.258862;  
//   rmsDeltaArray[3] =  1.3293;

//   diffErrArray[4] =  0;    
//   rDeltaArray[4] =  -0.164512; 
//   rmsDeltaArray[4] = 1.43615;

//   diffErrArray[5] =  0.05;   
//   rDeltaArray[5] =  -0.0674086;  
//   rmsDeltaArray[5] = 1.63776;

//   diffErrArray[6] = 0.1 ;
//   rDeltaArray[6] =  0.0284195;  
//   rmsDeltaArray[6] = 1.90095;

//   diffErrArray[7] = 0.15;   
//   rDeltaArray[7] =  0.122786;   
//   rmsDeltaArray[7] =  2.19241;

//   diffErrArray[8] =  0.2;
//    rDeltaArray[8] =  0.211219;
//    rmsDeltaArray[8] =  2.51358;
  

//   TCanvas *canHead = new TCanvas("canHead","canHead",600,400);
//   TGraph *rGraph = new TGraph(9,diffErrArray,rmsDeltaArray);
//   rGraph->Draw("ap");
//   rGraph->GetXaxis()->SetTitle("delta R (m)");
//   //  rGraph->GetYaxis()->SetTitle("Sum of (actual - expected) time (s)");
//   rGraph->GetYaxis()->SetTitle("Average RMS");
//   rGraph->SetTitle("Calculation of the sum of actual - expected time for all pairs");
//   rGraph->SetMarkerStyle(22);
//   rGraph->SetMarkerColor(1);
//   TF1 *f1 = new TF1("f1","x*[0]+[1]",-0.2,0.2);
//   f1->SetLineWidth(0.1);  
//   //rGraph->Fit("f1","R");
//   canHead->Update();
    
 //*********MINUIT METHOD*******************
   myMin->SetPrintLevel(-1);
   myMin->Migrad();   
   //myMin->GetParameter(0,deltaR,deltaRErr);
   //myMin->GetParameter(1,relDelta,relDeltaErr);
   //myMin->GetParameter(3,deltaZ,relDeltaZErr);
   //myMin->GetParameter(4,relDeltaZ,relDeltaZErr);
   //myMin->GetParameter(5,deltaPhi,relDeltaPhiErr);
   //myMin->GetParameter(6,relDeltaPhi,relDeltaPhiErr);
   
   for(int r = 0; r < 16; r++){
     myMin->GetParameter(0,deltaR[r],deltaRErr[r]);
   }
   
   for(int z = 0; z < 16; z++){
     // myMin->GetParameter(z+16,deltaZ[z],deltaZErr[z]);
   }

   for(int phi = 0; phi < 16; phi++){
      myMin->GetParameter(phi+1,deltaPhi[phi],deltaRErr[phi]);
   }

      myMin->GetParameter(18,deltaHeading,deltaHeading);

   //Double_t testIn1[7] = {deltaR,relDelta,antIn1,deltaZ,relDeltaZ,deltaPhi,relDeltaPhi};
   //Double_t testIn2[7] = {0,0,antIn1,0,0,0,0};

   //Double_t initialValue = antOpt(testIn2);
   //Double_t finalValue = antOpt(testIn1);

    //Double_t initialValue = 1;
    //Double_t finalValue = 1;

   //deltaROut = deltaR;
   //relDeltaOut = relDelta;

   //std::cout << "  " <<std::endl;
   //std::cout << "deltaR =  "  << deltaR << ", error =  " <<  deltaRErr  << ", initial value =  " << initialValue << ", final Value = " << finalValue << std::endl;
   //std::cout << "relDelta =  "  << relDelta << ", error =  " <<  relDeltaErr  << ", initial value =  " << initialValue << ", final Value = " << finalValue << std::endl;

   myMin->DeleteArrays();

   for(int r = 0; r < 16; r++){
     cout << deltaR[r] <<"  "<< 0 << "  "<< deltaPhi[r] << endl;
     optFileOut << deltaR[r] <<"  "<< 0 << "  "<< deltaPhi[r] << endl;
     cout << deltaHeading << endl;
   }   
   optFileOut.close();
   myMin->DeleteArrays();

  }

}

void iHateRoot(Int_t& npar, Double_t* gin,
                       Double_t& f, Double_t* par, Int_t flag){


  Double_t diffErr = antOpt(par);
  //Double_t diffErr = antOpt(par,headTree,adu5PatTree,corTree);
  f=diffErr;

}


