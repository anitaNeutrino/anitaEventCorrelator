
#ifndef FFTTOOLS_H
#define FFTTOOLS_H

//#include <fftw3.h>
#include "TGraph.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

//My includes
#include "FFTWComplex.h"


class FFTtools
{
public:
    FFTtools();
    ~FFTtools();
    
    static int runUnixCommandAndRedicrectOutputToTempFile(TString &theCommand, TString &theOutputFile);
    
    static double getAbs(FFTWComplex &theNum);
    static double *doInvFFT(int length, FFTWComplex *theInput);
    static FFTWComplex *doFFT(int length,double *theInput);
    
    static TGraph *makePowerSpectrum(TGraph *grWave);
    static TGraph *makeRawPowerSpectrum(TGraph *grWave);
    static TGraph *getCorrelationGraph(TGraph *gr1, TGraph *gr2);
    static double *getCorrelation(TGraph *gr1, TGraph *gr2,int firstIndex,int lastIndex);
    static double *getCorrelation(int length,float *oldY1, float *oldY2);
    static double *getCorrelation(int length,double *oldY1, double *oldY2);
    
    static TGraph *makeInverseInverseSpectrum(TGraph *grWave);
    static TGraph *combineGraphsUsingFFTs(Int_t numGraphs, TGraph **grPtr,double *theWeights=0);
    static Double_t *combineValuesUsingFFTs(Int_t numArrays, Double_t **thePtrPtr, Int_t eachLength);


};
   
#endif //FFTTOOLS_H
