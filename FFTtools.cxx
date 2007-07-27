#include "FFTtools.h"
#include <iostream>
#include <fstream>
#include "TSystem.h"

#include <fftw3.h>

using namespace std;

FFTtools::FFTtools()
{
}

FFTtools::~FFTtools()
{
}



FFTWComplex *FFTtools::doFFT(int length, double *theInput) {

//     TString theCommand("sillyFFT -f ");
//     char theNumber[40];
//     for(int i=0;i<length;i++) {
// 	sprintf(theNumber,"%f ",theInput[i]);
// //	cout << theNumber << "\t" << theInput[i] << endl;
// 	theCommand+=theNumber;
//     }

//     TString theOutputFile;
//     runUnixCommandAndRedicrectOutputToTempFile(theCommand,theOutputFile);
    
//     ifstream TheFile(theOutputFile.Data());
//     if(!TheFile) {
// 	return 0;
//     }
       
//     FFTWComplex *theOutput = new FFTWComplex [(length/2)+1];
//     int uptoNum=0;
//     double real,imag;
//     while((TheFile >> real >> imag) ) {
// 	if(uptoNum<(length/2)+1) {
// 	    theOutput[uptoNum].re=real;
// 	    theOutput[uptoNum].im=imag;
// 	    uptoNum++;
// 	}
// 	else break;
//     }
//     if(uptoNum==(length/2)+1) 
// 	return theOutput;
//     return 0;



//Here is what the sillyFFT program should be doing;    
    fftw_complex *theOutput = new fftw_complex [(length/2)+1];
    double *newInput = new double [length]; 
    
//(FFTWComplex*) fftw_malloc(sizeof(FFTWComplex)*((length/2)+1));
    //cout << length  << " " << input[0] << " " << theFFT[0] << endl;
    fftw_plan thePlan = fftw_plan_dft_r2c_1d(length,newInput,theOutput,FFTW_MEASURE);
    if(!thePlan) {
	cout << "Bollocks" << endl;
    }
        
    for(int i=0;i<length;i++) {
	newInput[i]=theInput[i];
    }

    for (int i=0;i<(length/2)+1;i++) {
	theOutput[i][0]=0.;
	theOutput[i][1]=0.;
    }
    fftw_execute(thePlan);
    delete [] newInput; 
    fftw_destroy_plan(thePlan);
    

    FFTWComplex *myOutput = new FFTWComplex [(length/2)+1];
    for (int i=0;i<(length/2)+1;i++) {
	myOutput[i].re=theOutput[i][0];
	myOutput[i].im=theOutput[i][1];
    }
    delete [] theOutput;
    return myOutput;
}


double *FFTtools::doInvFFT(int length, FFTWComplex *theInput) {
//    cout << "here in doInvFFT" << endl;

//     Int_t complexLength=(length/2)+1;

//     TString theCommand("sillyFFT -i ");
//     char theNumber[40];
//     for(int i=0;i<complexLength;i++) {
// 	sprintf(theNumber,"%e ",theInput[i].re);
// 	theCommand+=theNumber;
// 	sprintf(theNumber,"%e ",theInput[i].im);
// //	cout << theInput[i].im << endl;
// //	cout << theNumber << endl;
// 	theCommand+=theNumber;
//     }
// //    cout << "done something" << endl;
// //    cout << theCommand.Data() << endl;
//     TString theOutputFile;
//     runUnixCommandAndRedicrectOutputToTempFile(theCommand,theOutputFile);
    
// //    cout << "doneCommand " << endl;
//     ifstream TheFile(theOutputFile.Data());
//     if(!TheFile) {
// 	return 0;
//     }
       
//     double *theOutput = new double [length];
//     int uptoNum=0;
//     double theNum;
//     while(TheFile >> theNum) {
// //	cout << theNum << endl;
// 	if(uptoNum<length) {
// 	    theOutput[uptoNum]=theNum;
// 	    uptoNum++;
// 	}
// 	else break;
//     }
//     if(uptoNum==length) 
// 	return theOutput;
//     return 0;

// This is what sillyFFT should be doing
 //    //Takes account of normailisation 
    fftw_complex *newInput = new fftw_complex [(length/2)+1];
    double *theOutput = new double [length]; 
    fftw_plan thePlan = fftw_plan_dft_c2r_1d(length,newInput,theOutput,FFTW_MEASURE);
           
    for(int i=0;i<((length/2)+1);i++) {
	newInput[i][0]=theInput[i].re;
	newInput[i][1]=theInput[i].im;
    }
    for(int i=0;i<length;i++) {
	theOutput[i]=0;
    }
    
    fftw_execute(thePlan);
    delete [] newInput; 
    fftw_destroy_plan(thePlan);
    for(int i=0;i<length;i++) {
	theOutput[i]/=length;
    }
    return theOutput;
}




Double_t *FFTtools::combineValuesUsingFFTs(Int_t numArrays, Double_t **thePtrPtr, Int_t eachLength) {
    FFTWComplex **theFFTs = new FFTWComplex* [numArrays];
    for(int i=0;i<numArrays;i++) {
	theFFTs[i]=doFFT(eachLength,thePtrPtr[i]);
    }

    int fftLength=(eachLength/2)+1;
    FFTWComplex *combinedFFT = new FFTWComplex [fftLength];


    for(int i=0;i<fftLength;i++) {
	double tempAbs0=getAbs(theFFTs[0][i]);
	double tempTotAbs=tempAbs0;
	for(int arNum=1;arNum<numArrays;arNum++) {
	    tempTotAbs+=getAbs(theFFTs[arNum][i]);
	}

	combinedFFT[i].re=theFFTs[0][i].re*(tempTotAbs/(tempAbs0*double(numArrays)));
	combinedFFT[i].im=theFFTs[0][i].im*(tempTotAbs/(tempAbs0*double(numArrays)));
    }

    for(int i=0;i<numArrays;i++) {
	delete [] theFFTs[i];
    }
    delete [] theFFTs;
    double *newValues=doInvFFT(eachLength,combinedFFT);
    delete [] combinedFFT;
    return newValues;
    
}


TGraph *FFTtools::combineGraphsUsingFFTs(Int_t numGraphs, TGraph **grPtr,double *theWeights) {

    double totalWeight=0;
    if(theWeights) {
	for(int i=0;i<numGraphs;i++) {
	    totalWeight+=theWeights[i];
//	    cout << "Weight " << i << "\t" << theWeights[i] << endl;
	}
    }

    FFTWComplex **theFFTs = new FFTWComplex* [numGraphs];
    for(int i=0;i<numGraphs;i++) {
	double *oldY=grPtr[i]->GetY();
	int oldLength=grPtr[i]->GetN();
	theFFTs[i]=doFFT(oldLength,oldY);
    }

    int fftLength=((grPtr[0]->GetN())/2)+1;
    FFTWComplex *combinedFFT = new FFTWComplex [fftLength];

    for(int i=0;i<fftLength;i++) {
	if(theWeights) {
	    double tempAbs0=getAbs(theFFTs[0][i]);
	    double tempTotAbs=tempAbs0*theWeights[0];
	    for(int grNum=1;grNum<numGraphs;grNum++) {
		tempTotAbs+=getAbs(theFFTs[grNum][i])*theWeights[grNum];
	    }
	    
	    combinedFFT[i].re=theFFTs[0][i].re*(tempTotAbs/(tempAbs0*double(totalWeight)));
	    combinedFFT[i].im=theFFTs[0][i].im*(tempTotAbs/(tempAbs0*double(totalWeight)));
	}
	else {
	    double tempAbs0=getAbs(theFFTs[0][i]);
	    double tempTotAbs=tempAbs0;
	    for(int grNum=1;grNum<numGraphs;grNum++) {
		tempTotAbs+=getAbs(theFFTs[grNum][i]);
	    }
	    
	    combinedFFT[i].re=theFFTs[0][i].re*(tempTotAbs/(tempAbs0*double(numGraphs)));
	    combinedFFT[i].im=theFFTs[0][i].im*(tempTotAbs/(tempAbs0*double(numGraphs)));
	}
    }

    for(int i=0;i<numGraphs;i++) {
	delete [] theFFTs[i];
    }
    delete [] theFFTs;

    double *newX=grPtr[0]->GetX();
    int newLength=grPtr[0]->GetN();
    double *newY=doInvFFT(newLength,combinedFFT);
    TGraph *grOut = new TGraph(newLength,newX,newY);
    delete [] combinedFFT;
    return grOut;
    
}



TGraph *FFTtools::getCorrelationGraph(TGraph *gr1, TGraph *gr2) {

    int length=gr1->GetN();
    int length2=gr2->GetN();

    int N=length;
    if(length<length2)
       N=length2;

    double *oldY1 = new double [N];
    double *oldY2 = new double [N];
    
    double x,y;
    for(int i=0;i<N;i++) {
       y=0;
       if(i<length) {
	  gr1->GetPoint(i,x,y);
       }
       oldY1[i]=y;
       y=0;
       if(i<length2) {
	  gr2->GetPoint(i,x,y);
       }
       oldY2[i]=y;
       
    }
    gr1->GetPoint(1,x,y);
    double deltaT=x;
    gr1->GetPoint(N/2,x,y);
    double offset=x;
    gr1->GetPoint(0,x,y);
    deltaT-=x;
    offset-=x;
    offset*=-1;
    double waveOffset=x;
    gr2->GetPoint(0,x,y);
    waveOffset-=x;
    offset+=waveOffset;

    double *xVals = new double [N];
    double *yVals = new double [N];
    double *corVals=getCorrelation(N,oldY1,oldY2);
    for(int i=0;i<N;i++) {
       xVals[i]=offset+i*deltaT;
       if(i<N/2) {
	  yVals[i]=corVals[i+N/2];
       }
       else {
	  yVals[i]=corVals[i-N/2];
       }
    }


    TGraph *grCor = new TGraph(N,xVals,yVals);
    delete [] oldY1;
    delete [] oldY2;
    delete [] xVals;
    delete [] yVals;
    delete [] corVals;
    
    return grCor;
}


double *FFTtools::getCorrelation(int length,float *oldY1, float *oldY2) 
{
    double *newY1 = new double [length];
    double *newY2 = new double [length];
    for(int i=0;i<length;i++) {
	newY1[i]=(double)oldY1[i];
	newY2[i]=(double)oldY2[i];
    }

    double *theCorr=getCorrelation(length,newY1,newY2);
    delete [] newY1;
    delete [] newY2;
    return theCorr;
}


double *FFTtools::getCorrelation(int length,double *oldY1, double *oldY2) 
{

//    cout << "Here in getCorrelation" << endl;
    FFTWComplex *theFFT1=doFFT(length,oldY1);
    FFTWComplex *theFFT2=doFFT(length,oldY2);
    

    int newLength=(length/2)+1;
//    cout << "newLength " << newLength << endl;
    FFTWComplex *tempStep = new FFTWComplex [newLength];
    int no2=length>>1;
    for(int i=0;i<newLength;i++) {
	double reFFT1=theFFT1[i].re;
	double imFFT1=theFFT1[i].im;
	double reFFT2=theFFT2[i].re;
	double imFFT2=theFFT2[i].im;

	//Real part of output 
	tempStep[i].re=(reFFT1*reFFT2+imFFT1*imFFT2)/double(no2);
	//Imaginary part of output 
	tempStep[i].im=(imFFT1*reFFT2-reFFT1*imFFT2)/double(no2);
    }
//    cout << "finished messing around" << endl;
    double *theOutput=doInvFFT(length,tempStep);
//    cout << "got inverse" << endl;
    delete [] theFFT1;
    delete [] theFFT2;
    delete [] tempStep;
    return theOutput;

}


double *FFTtools::getCorrelation(TGraph *gr1, TGraph *gr2,int firstIndex,int lastIndex) {
    int tempLength=gr1->GetN();
    if(firstIndex<0 || lastIndex>tempLength) return 0;
    
    int length=lastIndex-firstIndex;
//    double *x1 = gr1->GetX();
    double *y1 = gr1->GetY();
//    double *x2 = gr2->GetX();
    double *y2 = gr2->GetY();
//    TGraph newGr1(length,&x1[firstIndex],&y1[firstIndex]);
//    TGraph newGr2(length,&x2[firstIndex],&y2[firstIndex]);
    
//    return getCorrelation(&newGr1,&newGr2);
    return getCorrelation(length,&y1[firstIndex],&y2[firstIndex]);
}


TGraph *FFTtools::makeInverseInverseSpectrum(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
//    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);
    double *invInvSpectrum = doInvFFT(length,theFFT);

    TGraph *grInvInv = new TGraph(length,oldX,invInvSpectrum);
//     for(int i=0;i<length;i++) {
// 	cout << oldX[i] << "\t" << invInvSpectrum[i] << endl;
//     }
    return grInvInv;

}



TGraph *FFTtools::makePowerSpectrum(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;
//    float deltaF=1.0/float(2*(newLength+2));
    double *newY = new double [newLength];
    double *newX = new double [newLength];

    double fMax = 1/(2*deltaT);  // In GHz
    double deltaF=fMax/newLength;
//    cout << "F max: " <<fMax << endl;
    //cout << deltaF << "\t" << otherDeltaF << endl;

    double tempF=0;
    double dt2df=deltaT/(newLength+2);// This is dt^2*df we need to normalize FFT power
    for(int i=0;i<newLength;i++) {
 	float power=pow(getAbs(theFFT[i]),2)*dt2df;
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power/=1e9; //convert to V^2/MHz
	if (power>0 ) power=10*TMath::Log10(power);
	else power=-100; //no reason
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }



    TGraph *grPower = new TGraph(newLength,newX,newY);
    return grPower;

}


TGraph *FFTtools::makeRawPowerSpectrum(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;
//    float deltaF=1.0/float(2*(newLength+2));
    double *newY = new double [newLength];
    double *newX = new double [newLength];

    double fMax = 1/(2*deltaT);  // In GHz
    double deltaF=fMax/newLength;
//    cout << "F max: " <<fMax << endl;
    //cout << deltaF << "\t" << otherDeltaF << endl;

    double tempF=0;
    double dt2df=deltaT/(newLength+2);// This is dt^2*df we need to normalize FFT power
    for(int i=0;i<newLength;i++) {
 	float power=pow(getAbs(theFFT[i]),2)*dt2df;
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
//	power/=1e9; //convert to V^2/MHz
	//if (power>0 ) power=10*TMath::Log10(power);
	//else power=-100; //no reason
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }



    TGraph *grPower = new TGraph(newLength,newX,newY);
    return grPower;

}


double FFTtools::getAbs(FFTWComplex &theNum) {
    return sqrt(theNum.re*theNum.re+theNum.im*theNum.im);
}


int FFTtools::runUnixCommandAndRedicrectOutputToTempFile(TString &theCommand, TString &theOutputFile ) {
    
    char *theFilename=gSystem->ConcatFileName(gSystem->TempDirectory(),
					      "bastardXXXXXX");
    int didItWork=mkstemp(theFilename);
    //cout << theFilename << endl;
    if(didItWork) {
	theOutputFile=theFilename;
	theCommand+=TString(" > ")+theFilename;
//	cout << theCommand << endl;
	delete theFilename;
	return gSystem->Exec(theCommand.Data());
    }
    else {
	delete theFilename;
	return 1;
    }


}
