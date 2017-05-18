//////////////////////////////////////////////////////////////////////////////
/////  PrettyAnitaEvent.h        Useful ANITA event class                      /////
/////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for plotting stuff like waveforms and correlations/////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////


/*! \mainpage ANITA Event Correlator
 *
 * \section intro_sec Introduction
 *
 * This is the somewhat sketchy documentation for my poor attempt at an event correlator library. At the moment this library is in pale in comparison to the work done by Andres, Jiwoo, Stephen and just about every other ANITA collaborator I suspect. But I did it so I'm sticking it here.
 *
 * \section prereq_sec Prerequisites
 *
 *  -# <A HREF="http://root.cern.ch">ROOT</A>
 *  -# <A HREF="http://www.fftw.org/">FFTW 3 -- Fastest Fourier Transform in the West</a>
 *  -# <A HREF="http://www.hep.ucl.ac.uk/uhen/anita/libRootFftwWrapper">libRootFftwWrapper -- a ROOT wrapper for FFTW 3</a>
 *  -# <A HREF="http://www.hep.ucl.ac.uk/uhen/anita/eventReader">ANITA-II Event Reader</a>
 * 
 * \section install_sec Installation
 * -# Checkout the code from the SVN repository, eg.: <BR><PRE>svn co https://delos.mps.ohio-state.edu/anitaGround/eventCorrelator/trunk myEventCorrelatorDir</PRE>
 * -# Define the ANITA_UTIL_INSTALL_DIR to point to the location you want the library installed (the library files will end up in (ANITA_UTIL_INSTALL_DIR)/lib and the header files in (ANITA_UTIL_INSTALL_DIR)/include).
 * -# Do <PRE>make</PRE><PRE>make install</PRE>
 * \section manual_sec Manual
 * If you are averse to reading web pages (and who wouldn't be) you can download a <a href="manual/eventCorrelator.pdf">pdf copy of the reference material</a> but be warned it won't be a thrilling read.
 */


#ifndef PRETTYANITAEVENT_H
#define PRETTYANITAEVENT_H
#include "UsefulAnitaEvent.h"
#include "CorrelationSummary.h"
#include "CorrelationSummaryAnita3.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <string>

#include <iostream>
#include <fstream>


//ANITA Includes
#include "CalibratedAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "FFTtools.h"
#include "FFTWComplex.h"


//ROOT Includes
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TVirtualFFT.h"
#include "TMinuit.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"





class CalibratedAnitaEvent;

//!  This is the event class, inherited from UsefulAnitaEvent, that has a number of correlation related methods.
/*!
  Basically this class is a good chunk of what the event correlator is all about. It inherits from UsefulAnitaEvent and can be used as a drop in replacement for the former. It provides a number of methods that do useful correlator-y stuff.
*/
class PrettyAnitaEvent: public UsefulAnitaEvent
{


 public:

  //! The assignment constructor.
  /*!
    \param eventPtr A pointer to the RawAnitaEvent.
    \param calType The desired calibration option.
    \param theHk The hk (needed for the temperature correction of the SURF timebase).
  */
   PrettyAnitaEvent(CalibratedAnitaEvent *calPtr, WaveCalType::WaveCalType_t calType=WaveCalType::kDefault);
   PrettyAnitaEvent(RawAnitaEvent *eventPtr,WaveCalType::WaveCalType_t calType, PrettyAnitaHk *theHk);
   
   PrettyAnitaEvent(RawAnitaEvent *eventPtr,WaveCalType::WaveCalType_t calType, RawAnitaHeader *headPtr);


  //Putative Analysis methods
  //! Calls getMaxAntennaCorrelation
  /*!
    \param pol Which polarisation to use?
    \param peakPtr An optional pointer to a double in which to store the peak/rms value.
    \return The antenna number of the antenna with the largest signal.
  */
  int getMaxAntenna(AnitaPol::AnitaPol_t pol, Double_t *peakPtr=0);
  //! Select the antenna with the maximum voltage squared.
  /*!
    \param pol Which polarisation to use?
    \param peakPtr An optional pointer to a double in which to store the peak V^2 value.
    \return The antenna number of the antenna with the largest V^2.
  */
  int getMaxAntennaVSquared(AnitaPol::AnitaPol_t pol, Double_t *peakPtr=0);
  //! Select the upper antenna with the maximum correlation (defined as peak/rms of the correlation) with it's pair in the lower ring.
  /*!
    \param pol Which polarisation to use?
    \param peakPtr An optional pointer to a double in which to store the peak/rms correlation value.
    \return The antenna number of the antenna with the largest peak/rms correlation value.
  */  
  int getMaxAntennaCorrelation(AnitaPol::AnitaPol_t pol, Double_t *peakPtr=0);
  //! Select the upper antenna with the maximum correlation (defined as peak/rms of the correlation) with it's pair in the lower ring.
  /*!
    \param pol Which polarisation to use?
    \param peakPtr An optional pointer to a double in which to store the peak/rms correlation value.
    \return The antenna number of the antenna with the largest peak/rms correlation value.
  */  
  int getMaxAntennaCorrelationRollingAvg(AnitaPol::AnitaPol_t pol, Double_t *peakPtr=0);

  //! Generates a CorrelationSummary object for a set of 10 antennas.
  /*!
    \param centreAnt The number of one of the antennas in the centre of the set of 10.
    \param pol Which polarisation to use?
    \param deltaT An optional value to use if interpolation is required. This value is taking as being the desired sampling period of the interpolated waveforms.
    \return A pointer to the CorrelationSummary object that is created.
  */  
  CorrelationSummary *getCorrelationSummary(Int_t centreAnt,AnitaPol::AnitaPol_t pol,Double_t deltaT=0);

  //! Generates a CorrelationSummaryAnita3 object for a set of 15 antennas.
  /*!
    \param centreAnt The number of one of the antennas in the centre of the set of 10.
    \param pol Which polarisation to use?
    \param deltaT An optional value to use if interpolation is required. This value is taking as being the desired sampling period of the interpolated waveforms.
    \return A pointer to the CorrelationSummaryAnita3 object that is created.
  */  
  CorrelationSummaryAnita3 *getCorrelationSummaryAnita3(Int_t centreAnt,AnitaPol::AnitaPol_t pol,Double_t deltaT=0);

  //! Generates a CorrelationSummaryAnita3 object for a set of 15 antennas.
  /*!
    \param centreAnt The number of one of the antennas in the centre of the set of 10.
    \param pol Which polarisation to use?
    \param deltaT An optional value to use if interpolation is required. This value is taking as being the desired sampling period of the interpolated waveforms.
    \return A pointer to the CorrelationSummaryAnita3 object that is created.
    Only the antenna positions are filled in this function, not the correlation values!!!
  */
  CorrelationSummaryAnita3 *createCorrelationSummaryAnita3(Int_t centreAnt,AnitaPol::AnitaPol_t pol, Double_t deltaT=0);

  //Canvas panel getters
  //! Generates a TCanvas with six waveforms plotted in it.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 6.
    \param pol Which polarisation to use?
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the waveforms.
  */  
  TCanvas *getSixWaveformCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  //! Generates a TCanvas with ten waveforms plotted in it.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 10.
    \param pol Which polarisation to use?
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the waveforms.
  */  
  TCanvas *getTenWaveformCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);

  //! Generates a TCanvas with six FFT power spectral densitity plots.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 6.
    \param pol Which polarisation to use?
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the PSDs plotted.
  */   
  TCanvas *getSixFFTPowerCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);

  //! Generates a TCanvas with six power envelope plots using FFTtools::getSimplePowerEnvelopeGraph
  /*!
    \param ant The number of one of the antennas in the centre of the set of 6.
    \param pol Which polarisation to use?
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the power enevlope graphs.
  */  
  TCanvas *getSixPowerEnvelopeCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  
  //! Generates a TCanvas with six interpolated waveforms plotted in it.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 6.
    \param pol Which polarisation to use?
    \param deltaT The sampling period of the interpolated waveforms (interpolation done using FFTtools::getInterpolatedGraph)
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the interpolated waveforms.
  */  
  TCanvas *getSixInterpolatedCanvas(int ant, AnitaPol::AnitaPol_t pol, Double_t deltaT=(1./(2.6*8)), TCanvas *can=0);
//! Generates a TCanvas with six correlations plotted in it.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 6.
    \param pol Which polarisation to use?
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the interpolated waveforms.
  */  
  TCanvas *getSixCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  //! Generates a TCanvas with ten correlations plotted in it.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 10.
    \param pol Which polarisation to use?
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the interpolated waveforms.
  */  
  TCanvas *getTenCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  //! Generates a TCanvas with eleven correlations (the combinations of six antennas) plotted in it.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 6.
    \param pol Which polarisation to use?
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the interpolated waveforms.
  */  
  TCanvas *getElevenCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol, TCanvas *can=0);
  //! Generates a TCanvas with eleven correlations of interpolated waveforms (the 11 correlations are combinations of six antennas) plotted in it.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 6.
    \param pol Which polarisation to use?
    \param deltaT The optional sampling period for the interpolated waveforms (default is up sampling by afactor of 8).
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the interpolated waveforms.
  */  
  TCanvas *getElevenInterpolationCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol,Double_t deltaT=(1./(2.6*8)), TCanvas *can=0);
  //! Generates a TCanvas with six correlations of interpolated waveforms plotted in it.
  /*!
    \param ant The number of one of the antennas in the centre of the set of 6.
    \param pol Which polarisation to use?
    \param deltaT The optional sampling period for the interpolated waveforms (default is up sampling by afactor of 8).
    \param can An optional pointer to a TCanvas to use instead of creating a new one.
    \return A pointer to the TCanvas with the interpolated waveforms.
  */  
  TCanvas *getSixInterpolatedCorrelationCanvas(int ant, AnitaPol::AnitaPol_t pol,Double_t deltaT=(1./(2.6*8)), TCanvas *can=0);
  //TCanvas *getTwoAntMap(int ant1,int ant2);


  //Graph getters
  TGraph *getSimplePowerEnvelopeGraph(int chanIndex); ///< Wrapper around FFTtools::getSimplePowerEnvelopeGraph
  TGraph *getInterpolatedGraph(int chanIndex, double deltaT); ///< Wrapper around FFTtools::getInterpolatedGraph
  TGraph *getFFTMagnitude(TGraph *grIn); ///< Wrapper around FFTtools::makePowerSpectrum
  TGraph *getFFTMagnitude(int chanIndex); ///< Wrapper around FFTtools::makePowerSpectrum
  TGraph *getCorrelation(int chanIndex1, int chanIndex2); ///< Wrapper around FFTtools::getCorrelationGraph
  TGraph *getCorrelation(TGraph *gr1, TGraph *gr2); ///< Wrapper around FFTtools::getCorrelationGraph
  TGraph *getCorrelationInterpolated(int chanIndex1, int chanIndex2, Double_t deltaT=(1./(2.6*8)) ); ///< Wrapper around FFTtools::getInterpolatedCorrelationGraph
  
  void fillSixAntArrays(int ant, int topAnts[3], int bottomAnts[3]); ///< Utility to get neighbouring antenna numbers
  void fillNextFourAntArrays(int ant, int nextFourAnts[4]);///< Utility to get next to neighbouring antenna numbers
  void fillNadirArrays(int ant, int nadirAnts[9]);

  void fillNineAntArrays(int ant, int nineAnts[9]); ///< Utility to get neighbouring antenna numbers ( Top 0-2, Middle 3-5, Bottom 6-8)
  void fillNextSixAntArrays(int ant, int nextFourAnts[4]);///< Utility to get next to neighbouring antenna numbers


#ifndef PLEASE_IGNORE_ME_DOXYGEN
  ClassDef(PrettyAnitaEvent,1); ///< ROOT's magic macro
#endif

  void setPassBandFilterFlag( int flag) { fPassBandFilter=flag;}
  void setNotchFilterFlag( int numNotches) { fNotchFilter=numNotches;}
  void setPassBandLimits(Double_t low, Double_t high)
     { fLowPassEdge=low; fHighPassEdge=high;}
  void setNotchBandLimits(Int_t notchNum, Double_t low, Double_t high)
     { fLowNotchEdge[notchNum]=low; fHighNotchEdge[notchNum]=high;}


 private:
  void setStyleSixCanvas(); ///< gStyle setup
  int getPrettyColour(int index); ///< Utility to get a pretty colour.

  Double_t fDeltaT; ///< The interpolated sampling rate.
  Double_t fWaveOffset; ///< The difference in T0 of two channels.
  Int_t fPassBandFilter; ///< Whether or not to pass band filter the interpolated waves;
  Int_t fNotchFilter; ///< Whether or not to notch filter, and how many notches
  Double_t fLowPassEdge; ///< The lower edge of the pass band
  Double_t fHighPassEdge; ///< The higher edge of the pass band
  Double_t fLowNotchEdge[10]; ///< The lower edge of the notch band
  Double_t fHighNotchEdge[10]; ///< The higher edge of the notch band
  

};


#endif //PRETTYANITAEVENT_H
