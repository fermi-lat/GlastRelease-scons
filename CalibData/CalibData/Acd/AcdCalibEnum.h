#ifndef CalibData_AcdCalibEnum_h
#define CalibData_AcdCalibEnum_h 

/**
 * @brief Data for ACD calibrations
 *
 * @author Eric Charles
 * $Header$
 **/

namespace AcdCalibData {

  /// Types of calibrations
  enum CALTYPE{NONE=-1,
	       PEDESTAL=0,        // Pedestal, stored in DB
	       GAIN=1,            // Gains (aka mips peaks), stored in DB
	       VETO=2,            // Veto thresholds, stored in DB
	       RANGE=3,           // Range crossover, stored in DB
	       CNO=4,	          // CNO thresholds, stored in DB
	       HIGH_RANGE=5,      // High Range calibration, stored in DB
	       COHERENT_NOISE=6,  // Coherent Noise calibration, stored in DB
	       TIME_PROF=7,       // Time profile plots, used in EMI testing
	       UNPAIRED=8,        // Looking for channels w/ only 1 PMT
	       HITMAP=9,          // Checking hitmap timing and latching
	       MERITCALIB=10,     // Filling an ntuple with variables from merit file
	       NDESC=11};  
 
};

#endif
