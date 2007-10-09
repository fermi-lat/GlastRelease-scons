#ifndef CalibData_AcdCalibEnum_h
#define CalibData_AcdCalibEnum_h 

namespace AcdCalibData {

  // Types of calibrations
  enum CALTYPE{NONE=-1,
	       PEDESTAL=0, 
	       GAIN=1, 
	       VETO=2,
	       RANGE=3,
	       CNO=4,	   
	       HIGH_RANGE=5,
	       COHERENT_NOISE=6,
	       TIME_PROF=7,
	       UNPAIRED=8, 
	       HITMAP=9,
	       MERITCALIB=10,
	       NDESC=11};  
 
};

#endif
