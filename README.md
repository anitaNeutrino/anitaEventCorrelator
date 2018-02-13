# anitaEventCorrelator

//////////////////////////////////////////////////////////////////////////////
////    eventCorrelator   -- December 2007                                ////
////   Ryan Nichol (rjn@hep.ucl.ac.uk)                                   /////
//////////////////////////////////////////////////////////////////////////////

Introduction
==============

This is the new eventCorrelator it relies on libRootFftWrapper to do the correlations. The way that it expects to work is having ANITA_UTIL_INSTALL_DIR set and it expects to find libRootFftWraper.so in $(ANITA_UTIL_INSTALL_DIR)/lib. Of course any location that LD and gcc looks will work. 


NOTE - for RampdemReader to work, ramp1kmdem_wgs_v2.hdr, ramp1kmdem_wgs_v2.bin and surface_elevation.asc are required, get them from www.hep.ucl.ac.uk/~mottram/ and place them in your eventCorrelator/data directory
