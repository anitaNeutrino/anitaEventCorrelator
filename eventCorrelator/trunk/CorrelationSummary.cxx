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

CorrelationSummary::CorrelationSummary( )
{

}

CorrelationSummary::~CorrelationSummary( )
{

}
CorrelationSummary::CorrelationSummary( int teventNumber, int tcentreAnt, int tsixAnts[], double dt)
   : eventNumber(teventNumber),centreAntenna(tcentreAnt),deltaT(dt)
{
  for(int i=0;i<6;i++) 
    sixAnts[i]=tsixAnts[i];

}


void CorrelationSummary::fillErrorsAndFit()
{
   //At some stage this should probably do something, otherwise it won't be much cop.


}
