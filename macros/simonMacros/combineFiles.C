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

void combineFiles() {

 ifstream oldfile("oldSimonNumbers.txt");
 ifstream newfile("newSimonNumbers.txt");     

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
     


       while (!oldfile.eof()){
	 oldfile >> number;
	oldfile >> dt;
       	oldfile >> dr;
       	oldfile >> dphi;
	oldfile >> dz;
 	
	newfile >> number;
	cout << number << "  ";
	newfile >> dt2;
	cout << dt+dt2 << "  ";
	newfile >> dr2;
	cout << dr+dr2 << "  ";
	newfile >> dphi2;
	cout << dphi+dphi2 << "  ";
	newfile >> dz2;
 	cout << dz+dz2 << "  " << endl;



       }



  }

}
