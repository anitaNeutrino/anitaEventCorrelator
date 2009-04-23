#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <stdio.h>

 #include "TGraph.h"
 #include "TCanvas.h"
 #include "TFile.h"
 #include "TH1F.h"
 #include "TVirtualFitter.h"
 #include "TStyle.h"
 #include "TPaveText.h"

using namespace std;

void deltaPlotter() {

 ifstream oldfile("deltaNumbers.txt");

  if (oldfile.fail()) {
    cout << "oldSimonNumbers.txt" << "  file  NOT FOUND" << endl;
    return 1;
  }
  
  else{
    
    double dummy = 0;
    int number = 0;
    double dt = 0;
    double dr = 0;
    double dphi = 0;
    double dz = 0;
    double dt2 = 0;
    double dr2 = 0;
    double dphi2 = 0;
    double dz2 = 0;  

     cout << "opening file   " << "oldSimonNumbers.txt" <<endl;
     
     double array1[41] = {0};
     double array2[41] = {0};
     int i = 0;

       while (!oldfile.eof()){
	 oldfile >> number;
	oldfile >> dt;
       	oldfile >> dr;
       	oldfile >> dphi;
	oldfile >> dz;
       
	array1[i] = number;
	array2[i] = dz;
	i++;

       }


TGraph *tempAntGraph  = new TGraph(40,  array1, array2);

	
	tempAntGraph->Draw("ap");    

	tempAntGraph->SetMarkerStyle(27);
   
	


  }

}
