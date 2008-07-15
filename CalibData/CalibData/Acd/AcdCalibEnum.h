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
	       PEDESTAL=0,        // Pedestal, stored in DB, used in MOOT
	       GAIN=1,            // Gains (aka mips peaks), stored in DB, used in MOOT
	       VETO=2,            // Veto thresholds, stored in DB
	       RANGE=3,           // Range crossover, stored in DB
	       CNO=4,	          // CNO thresholds, stored in DB
	       HIGH_RANGE=5,      // High Range calibration, stored in DB
	       COHERENT_NOISE=6,  // Coherent Noise calibration, stored in DB
	       RIBBON=7,          // Ribbon light attenuation, stored in DB
	       PED_HIGH=8,        // High range pedestals, used in MOOT
	       CARBON=9,          // Carbon Peak, used in MOOT
	       VETO_FIT=10,       // Veto_dac -> PHA mapping, used in MOOT
	       CNO_FIT=11,        // hld_dac -> PHA mapping, used in MOOT
	       TIME_PROF=12,      // Time profile plots, used in EMI testing
	       UNPAIRED=13,       // Looking for channels w/ only 1 PMT
	       HITMAP=14,         // Checking hitmap timing and latching
	       MERITCALIB=15,     // Filling an ntuple with variables from merit file
	       NDESC=16};  
 
};

#endif
