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
  //Here is what the sillyFFT program should be doing;    
    fftw_complex *theOutput = new fftw_complex [(length/2)+1];
    double *newInput = new double [length]; 
    
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
   //Now we'll extend this up to a power of 2
    int length=gr1->GetN();
    int length2=gr2->GetN();

    int N=int(TMath::Power(2,int(TMath::Log2(length))+2));
    if(N<length2)
       N=int(TMath::Power(2,int(TMath::Log2(length2))+2));

    //Will really assume that N's are equal for now
    int firstRealSamp=(N-length)/2;

    double *oldY1 = new double [N];
    double *oldY2 = new double [N];
    
    double x,y;
    Double_t x2,y2;
    gr1->GetPoint(1,x2,y2);
    gr1->GetPoint(0,x,y);
    double deltaT=x2-x;
    double firstX=x;

    gr2->GetPoint(0,x2,y2);
    double waveOffset=firstX-x2;
    

    //    gr1->GetPoint(N/2,x2,y2);
    //    double offset=x-x2;
    //    std::cout << length << "\t" << length2 << "\n";

    for(int i=0;i<N;i++) {
       
       if(i<firstRealSamp || i>=firstRealSamp+length)
	  y=0;
       else {
	  gr1->GetPoint(i-firstRealSamp,x,y);
       }
       oldY1[i]=y;
	  
       if(i<firstRealSamp || i>=firstRealSamp+length2)
	  y=0;
       else {
	  gr2->GetPoint(i-firstRealSamp,x,y);
       }
       oldY2[i]=y;
              
    }


    //    offset+=waveOffset;

    double *xVals = new double [N];
    double *yVals = new double [N];
    double *corVals=getCorrelation(N,oldY1,oldY2);
    for(int i=0;i<N;i++) {
       if(i<N/2) {
	  //Positive
	  xVals[i+(N/2)]=(i*deltaT)+waveOffset;
	  yVals[i+(N/2)]=corVals[i];
       }
       else {
	  //Negative
	  xVals[i-(N/2)]=((i-N)*deltaT)+waveOffset;
	  yVals[i-(N/2)]=corVals[i];	  
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

TGraph *FFTtools::getHilbertTransform(TGraph *grWave)
{
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  //    double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();
  FFTWComplex *theFFT=doFFT(length,oldY);  
  int newLength=(length/2)+1;
  for(int i=0;i<newLength;i++) {
    double tempIm=theFFT[i].im;
    theFFT[i].im=theFFT[i].re;
    theFFT[i].re=-1*tempIm;
  }
  double *hilbert = doInvFFT(length,theFFT);
  
  TGraph *grHilbert = new TGraph(length,oldX,hilbert);
  delete [] hilbert;
  delete [] theFFT;

  return grHilbert;
}



TGraph *FFTtools::makePowerSpectrum(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;
    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In GHz
    double deltaF=1/(deltaT*length);

    double tempF=0;
    for(int i=0;i<newLength;i++) {
       float power=pow(getAbs(theFFT[i]),2);
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power/=length;
	//	if (power>0 ) power=10*TMath::Log10(power);
	//	else power=-100; //no reason
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }



    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] newY;
    delete [] newX;
    return grPower;

}



TGraph *FFTtools::makePowerSpectrumPeriodogram(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);
    
    int newLength=(length/2)+1;

    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length);

    double tempF=0;
    for(int i=0;i<newLength;i++) {
       float power=pow(getAbs(theFFT[i]),2);
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power/=double(length*length);
	//	if (power>0 ) power=10*TMath::Log10(power);
	//	else power=-100; //no reason
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }

    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] newY;
    delete [] newX;
    return grPower;

}

TGraph *FFTtools::makePowerSpectrumVoltsSeconds(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;

    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e-6; //MHz


    double tempF=0;
    for(int i=0;i<newLength;i++) {
       float power=pow(getAbs(theFFT[i]),2);
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power*=deltaT/(length); //For time-integral squared amplitude
	power/=deltaF;//Just to normalise bin-widths
	//Ends up the same as dt^2, need to integrate the power (multiply by df)
	//to get a meaningful number out.

	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }


    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] newY;
    delete [] newX;
    return grPower;

}


TGraph *FFTtools::makePowerSpectrumVoltsSecondsdB(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;

    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e-6; //MHz


    double tempF=0;
    for(int i=0;i<newLength;i++) {
       float power=pow(getAbs(theFFT[i]),2);
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power*=deltaT/(length); //For time-integral squared amplitude
	power/=deltaF;//Just to normalise bin-widths
	//Ends up the same as dt^2, need to integrate the power (multiply by df)
	//to get a meaningful number out.	
	
	if (power>0 ) power=10*TMath::Log10(power);
	else power=-1000; //no reason
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }


    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] newY;
    delete [] newX;
    return grPower;

}

TGraph *FFTtools::makePowerSpectrumVoltsSecondsPadded(TGraph *grWave, Int_t padFactor) {

   TGraph *grPad=padWave(grWave,padFactor);
   TGraph *grPower=makePowerSpectrumVoltsSeconds(grPad);
   delete grPad;
   return grPower;
   
}


TGraph *FFTtools::makePowerSpectrumVoltsSecondsPaddeddB(TGraph *grWave, Int_t padFactor) {
   TGraph *grPad=padWave(grWave,padFactor);
   TGraph *grPower=makePowerSpectrumVoltsSecondsdB(grPad);
   delete grPad;
   return grPower;
}


TGraph *FFTtools::makeRawPowerSpectrum(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;
    double *newY = new double [newLength];
    double *newX = new double [newLength];

    double deltaF=1/(deltaT*length);
    //    double fMax = 1/(2*deltaT);  // In GHz

    double tempF=0;
    for(int i=0;i<newLength;i++) {
       float power=pow(getAbs(theFFT[i]),2);
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }


    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] newY;
    delete [] newX;
    return grPower;

}


double FFTtools::getAbs(FFTWComplex &theNum) {
    return sqrt(theNum.re*theNum.re+theNum.im*theNum.im);
}


Double_t FFTtools::sumPower(TGraph *gr,Int_t firstBin,Int_t lastBin)
{
  Double_t integral=0;
  Double_t freq,power;
  if(firstBin<0) firstBin=0;
  if(lastBin<0) lastBin=gr->GetN()-1;
  for(int i=firstBin;i<=lastBin;i++) {
    gr->GetPoint(i,freq,power);
    integral+=power;
  }
  return integral;
}

Double_t FFTtools::integratePower(TGraph *gr,Int_t firstBin,Int_t lastBin)
{
  Double_t integral=0;
  Double_t freq,power;
  gr->GetPoint(1,freq,power);
  Double_t df=freq;
  gr->GetPoint(0,freq,power);
  df-=freq;
  if(firstBin<0) firstBin=0;
  if(lastBin<0) lastBin=gr->GetN()-1;
  for(int i=firstBin;i<=lastBin;i++) {
    gr->GetPoint(i,freq,power);
    integral+=power*df;
  }
  return integral;
}

Double_t FFTtools::sumVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin)
{
  Double_t integral=0;
  Double_t time,volts;
  if(firstBin<0) firstBin=0;
  if(lastBin<0) lastBin=gr->GetN()-1;
  for(int i=firstBin;i<=lastBin;i++) {
    gr->GetPoint(i,time,volts);
    integral+=volts*volts;
  }
  return integral;
}

Double_t FFTtools::integrateVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin)
{
  Double_t integral=0;
  Double_t time,volts;
  gr->GetPoint(1,time,volts);
  Double_t dt=time;
  gr->GetPoint(0,time,volts);
  dt-=time;
  if(firstBin<0) firstBin=0;
  if(lastBin<0) lastBin=gr->GetN()-1;
  for(int i=firstBin;i<=lastBin;i++) {
    gr->GetPoint(i,time,volts);
    integral+=volts*volts*dt;
  }
  return integral;
}

Int_t FFTtools::getPeakBin(TGraph *gr) 
{
  Double_t x,y;
  gr->GetPoint(0,x,y);
  Double_t peakVal=y;
  Int_t peakBin=0;
  for(int i=1;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    if(peakVal<y) {
      peakVal=y;
      peakBin=i;
    }      
  }
  return peakBin;
}

TGraph *FFTtools::subtractGraphs(TGraph *grA, TGraph *grB) 
{
  Int_t N1=grA->GetN();
  Int_t N2=grB->GetN();
  if(N1!=N2) return NULL;

  Double_t *newY = new Double_t [N1];
  Double_t *xVals=grA->GetX();
  Double_t x,yA,yB;
  for(int i=0;i<N1;i++) {
    grA->GetPoint(i,x,yA);
    grB->GetPoint(i,x,yB);
    newY[i]=yA-yB;
  }
  TGraph *grDiff = new TGraph(N1,xVals,newY);
  delete [] newY;
  return grDiff;
}

TGraph *FFTtools::divideGraphs(TGraph *grA, TGraph *grB) 
{
  Int_t N1=grA->GetN();
  Int_t N2=grB->GetN();
  if(N1!=N2) return NULL;

  Double_t *newY = new Double_t [N1];
  Double_t *xVals=grA->GetX();
  Double_t x,yA,yB;
  for(int i=0;i<N1;i++) {
    grA->GetPoint(i,x,yA);
    grB->GetPoint(i,x,yB);
    newY[i]=yA/yB;
  }
  TGraph *grRat = new TGraph(N1,xVals,newY);
  delete [] newY;
  return grRat;
}


TGraph *FFTtools::ratioSubtractOneGraphs(TGraph *grA, TGraph *grB) 
{
   Int_t N1=grA->GetN();
  Int_t N2=grB->GetN();
  //  if(N1!=N2) return NULL;
  
  Int_t newN=N1;
  if(N2<N1) {
     newN=N2;
     //return NULL;
  }
  Double_t *xVals=grA->GetX();
  Double_t *xBVals=grB->GetX();
  Double_t deltaF=xVals[1]-xVals[0];
  Double_t deltaFB=xBVals[1]-xBVals[0];
  
  if(TMath::Abs(deltaFB-deltaF)>1) return NULL;
  //  cout << newN << endl;
  Double_t *newY = new Double_t [newN];
  Double_t x,yA,yB;
  for(int i=0;i<newN;i++) {
    grA->GetPoint(i,x,yA);
    grB->GetPoint(i,x,yB);
    newY[i]=1-yA/yB;
  }
  TGraph *grRat = new TGraph(newN,xVals,newY);
  delete [] newY;
  return grRat;
}

TGraph *FFTtools::dbGraphs(TGraph *grA, TGraph *grB) 
{
  Int_t N1=grA->GetN();
  Int_t N2=grB->GetN();
  //  if(N1!=N2) return NULL;
  
  Int_t newN=N1;
  if(N2<N1) {
     newN=N2;
     //return NULL;
  }
  Double_t *xVals=grA->GetX();
  Double_t *xBVals=grB->GetX();
  Double_t deltaF=xVals[1]-xVals[0];
  Double_t deltaFB=xBVals[1]-xBVals[0];
  //  cout << N1 << "\t" << N2 << "\t" << deltaF << "\t" << deltaFB << "\n";


  if(TMath::Abs(deltaFB-deltaF)>1) return NULL;
  //  cout << newN << endl;
  Double_t *newY = new Double_t [newN];
  Double_t x,yA,yB;
  for(int i=0;i<newN;i++) {
    grA->GetPoint(i,x,yA);
    grB->GetPoint(i,x,yB);
    newY[i]=10*TMath::Log10(yA/yB);
  }
  TGraph *grRat = new TGraph(newN,xVals,newY);
  delete [] newY;
  return grRat;
}


TGraph *FFTtools::smoothFFT(TGraph *gr,Int_t factor) 
{
  Int_t N=gr->GetN();
  Int_t newN=N/factor;
  Double_t *xVals=gr->GetX();
  Double_t *yVals=gr->GetY();

  Double_t *newX = new Double_t [newN];
  Double_t *newY = new Double_t [newN];
  Double_t sumX=0;
  Double_t sumY=0;
  //  cerr << N << "\t" << factor << "\t" << newN << "\n";
  for(int i=0;i<N;i++) {
     sumX+=xVals[i];
     sumY+=yVals[i];
     if((i+1)%factor==0) {
	//	cerr << i << "\t" << sumX << "\t" << sumY << "\n";
	//New Point
	newX[(i+1)/factor-1]=sumX/factor;
	newY[(i+1)/factor-1]=sumY/factor;
	sumX=0;
	sumY=0;
     }
  }
  TGraph *grSmooth = new TGraph(newN,newX,newY);
//   delete [] newX;
//   delete [] newY;
  return grSmooth;
}

TGraph *FFTtools::padWave(TGraph *grWave, Int_t padFactor) {
   double *oldY = grWave->GetY();
   double *oldX = grWave->GetX();
   double deltaT=oldX[1]-oldX[0];
   int realLength = grWave->GetN();
   int length = grWave->GetN()*padFactor;
   double *paddedY = new double [length];
   double *paddedX = new double [length];
   int newStart=(realLength*(padFactor-1))/2;
   for(int i=0;i<length;i++) {
      int waveIndex=i-newStart;
      paddedY[i]=0;
      paddedX[i]=(waveIndex*deltaT)+oldX[0];
      if(waveIndex>=0 && waveIndex<realLength) {
	 paddedY[i]=oldY[waveIndex];
      }
   }
   TGraph *grPadded = new TGraph(length,paddedX,paddedY);
   delete [] paddedX;
   delete [] paddedY;
   return grPadded;
}

TGraph *FFTtools::getSimplePowerEnvelopeGraph(TGraph *gr) {
  Double_t *ySq = new Double_t [gr->GetN()];
  Double_t *xOrig = new Double_t [gr->GetN()];
  Double_t *yEnvelope = new Double_t[gr->GetN()];
  Double_t *xEnvelope = new Double_t[gr->GetN()];
  
  Double_t x,y;
  Int_t numPoints=0;
  
  for(int i=0;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    ySq[i]=y*y;
    xOrig[i]=x;
    if(i==1) {
      if(ySq[0]>ySq[i]) {
        yEnvelope[numPoints]=ySq[0];
        xEnvelope[numPoints]=xOrig[0];
        numPoints++;
      }
    }
    else if(i==gr->GetN()-1 && ySq[i]>ySq[i-1]) {
      yEnvelope[numPoints]=ySq[i];
      xEnvelope[numPoints]=xOrig[i];
      numPoints++;
    }
    else if(ySq[i-1]>ySq[i-2] && ySq[i-1]>ySq[i]) {
      yEnvelope[numPoints]=ySq[i-1];
      xEnvelope[numPoints]=xOrig[i-1];
      numPoints++;
    }
  }                                  
  TGraph *grEnvelope= new TGraph(numPoints,xEnvelope,yEnvelope);
  delete [] ySq;
  delete [] xOrig;
  delete [] yEnvelope;
  delete [] xEnvelope;
  return grEnvelope;
}
