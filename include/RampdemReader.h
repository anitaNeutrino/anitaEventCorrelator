/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt, Matt Mottram, Stephen Hoover, probably others.
 Email: strutt@physics.ucla.edu was the most recent contributor

 Description:
 Class to read in the RAMPDEM data... and now all the BEDMAP2 data sets.
***********************************************************************************************************/

#ifndef RAMPDEMREADER_H
#define RAMPDEMREADER_H

//Includes
#include <iostream>
#include "TObject.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TGaxis.h"
#include <vector>
#include <map>


/**
 * @class Class to read in the RAMPDEM data and now all the BEDMAP2 data sets in one place.
 *
 * Re-implemented from singleton members -> to static, file-only variables in the RampdemReader.cxx file.
 * Instance function can still create a single instance, but it's now pointless.
 * At some future point this class can be reimplemented as a namespace a la FFTtools.
 */
class RampdemReader{

public:


  // Enumerates the sets of rampdem/bedmap2 data
  // OK, these files are enormous! There's no way I'm committing them all to github
  // Originals are from here: https://secure.antarctica.ac.uk/data/bedmap2/
  // Feel free to download the rest if you want all the other fancy ones
  typedef enum{
    rampdem, // this is the old data in surfaceElevation.asc (a 6MB file), this is included by default.
    bed,
    // coverage,
    // grounded_bed_uncertainty,
    icemask_grounded_and_shelves,
    // lakemask_vostok,
    // rockmask,
    surface,
    thickness,
    // bedmap2_thickness_uncertainty_5km
  } dataSet;


  RampdemReader();
  ~RampdemReader();

  // Class could (should) now be implemented as a namespace with loads any data behind the scenes...
  // All the cool kids are marking things as deprecated these days
  static RampdemReader*  Instance() __attribute__((deprecated("All methods are now static, e.g. call with RampdemReader::someMethod()"))); ///<Instance generator

  static Double_t Geoid(Double_t latitude);
  static Double_t Area(Double_t latitude, RampdemReader::dataSet=rampdem);

  static void ENtoLonLat(Int_t e_coord, Int_t n_coord,Double_t& lon, Double_t& lat,
			 RampdemReader::dataSet=rampdem);

  static void LonLattoEN(Double_t lon, Double_t lat, int& e_coord, int& n_coord,
			 RampdemReader::dataSet=rampdem);

  static void EastingNorthingToEN(Double_t easting, Double_t northing, Int_t &e_coord, Int_t &n_coord,
				  RampdemReader::dataSet=rampdem);

  static void LonLatToEastingNorthing(Double_t lon, Double_t lat, Double_t &easting, Double_t &northing);


  // This one should be just a geometrical calculation... however, it uses goes via the EN bins of a dataset
  // and the internal representation changes between RampDem and BEDMAP2... so beware...
  // (Can't be bothered to change this at the present time)
  static void EastingNorthingToLonLat(Double_t easting, Double_t northing, Double_t &lon, Double_t &lat, RampdemReader::dataSet=rampdem);

  static Bool_t isOnContinent(Double_t lon, Double_t lat);

  //Data Output methods

  static Double_t Surface(Double_t longitude, Double_t latitude);
  static Double_t SurfaceAboveGeoid(Double_t longitude, Double_t latitude, RampdemReader::dataSet=rampdem);
  static Double_t SurfaceAboveGeoidRampDem(Double_t longitude, Double_t latitude);

  static TProfile2D *rampMap(int coarseness_factor, int set_log_scale,UInt_t &xBins,UInt_t &yBins) __attribute__((deprecated("Prefer RampdemReader::getMap() with RampdemReader::rampdem data set")));

  static TProfile2D *rampMapPartial(int coarseness_factor, double centralLon, double centralLat, double rangeMetres, Int_t &xBins, Int_t &yBins, Double_t &xMin, Double_t &xMax, Double_t &yMin, Double_t &yMax)__attribute__((deprecated("Prefer RampdemReader::getMapPartial() with RampdemReader::rampdem data set")));


  static TGaxis *distanceScale(Double_t xMin,Double_t xMax,Double_t yMin,Double_t yMax);

  //Generic method to flip Endianness.
  //WARNING: Flips byte order of anything put in - do not use on things like stuctures or classes!
  template <class thing>
  inline static void flipEndian(thing &in){
    int size = sizeof(thing);

    thing out;

    char* p_in = (char *) &in;
    char* p_out = (char *) &out;

    for(int i=0;i<size;i++)
      p_out[i] = p_in[size-1-i];

    in = out;

    return;
  }

  static void getMapCoordinates(double &xMin,double &yMin,double &xMax,double &yMax, RampdemReader::dataSet=rampdem);
  static TProfile2D* getMap(RampdemReader::dataSet dataSet, int coarseness_factor);
  static TProfile2D* getMapPartial(RampdemReader::dataSet dataSet, int coarseness, double centralLon, double centralLat, double rangeMetres);

  static int readRAMPDEM(); // no need to call this explicitly, is called once only when rampdem data is requested

  static const char* dataSetToAxisTitle(RampdemReader::dataSet dataSet);
  static void getNumXY(Int_t& numX, Int_t&numY, RampdemReader::dataSet dataSet=rampdem);

  static TProfile2D* fillThisHist(TProfile2D* theHist, RampdemReader::dataSet dataSet);
protected:
  static RampdemReader *fgInstance; // deprecated, see RampdemReader::Instance() function.
  // protect against multiple instances


};

#endif //RAMPDEMREADER_H
