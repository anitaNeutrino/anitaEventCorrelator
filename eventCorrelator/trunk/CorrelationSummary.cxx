//////////////////////////////////////////////////////////////////////////////
/////  CorrelationSummary.cxx        ANITA event reading class                  /////////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for plotting event stuff like waveforms and correlations/////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//ANITA Includes
#include "CorrelationSummary.h"
#include "AnitaGeomTool.h"

//ROOT Includes
#include "TMath.h"
#include "TStyle.h"



ClassImp(CorrelationSummary);

CorrelationSummary::CorrelationSummary( int teventNumber, int tmaxAnt, int tsixAnts[])
  : eventNumber(teventNumber),maxAntenna(tmaxAnt)
{
  for(int i=0;i<6;i++) 
    sixAnts[i]=tsixAnts[i];

}
