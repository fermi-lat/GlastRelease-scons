// $Header$

#ifndef CalibData_CalibModel_h
#define CalibData_CalibModel_h

/**
       @file CalibModel.h
 Definition of strings holding paths to calibration data objects
 in the TDDS (transient detector data store)
   @author : adapted from similar structures for TDS
   @author   J. Bogart
 */ 

#include <string>
#include <vector>
#include "GaudiKernel/ClassID.h"

// For everybody except the CalibModel class implementation file,
// the variables are extern.  CalibModel actually defines them.

#if defined(_CalibData_CalibModel_cxx)
#define  _EXTERN_ 
#else
#define  _EXTERN_ extern
#endif

    namespace CalibData {
      // ACD calib types
      _EXTERN_ std::string   ACD_Eff;
      _EXTERN_ std::string   ACD_ThreshHigh;
      _EXTERN_ std::string   ACD_ThreshVeto;
      _EXTERN_ std::string   ACD_Ped;
      _EXTERN_ std::string   ACD_ElecGain;

      // TKR calib types
      _EXTERN_ std::string   TKR_BadChan;
      _EXTERN_ std::string   TKR_HotChan;
      _EXTERN_ std::string   TKR_DeadChan;
      _EXTERN_ std::string   TKR_TOTSignal;
      _EXTERN_ std::string   TKR_TOTDist;
      _EXTERN_ std::string   TKR_MIPEff;

      _EXTERN_ std::string   CAL_LightAtt;
      _EXTERN_ std::string   CAL_LightAsym;
      _EXTERN_ std::string   CAL_LightYield;
      _EXTERN_ std::string   CAL_Ped;
      _EXTERN_ std::string   CAL_ElecGain;
      _EXTERN_ std::string   CAL_IntNonlin;
      _EXTERN_ std::string   CAL_DiffNonlin;
      _EXTERN_ std::string   CAL_HotChan;
      _EXTERN_ std::string   CAL_DeadChan;
      _EXTERN_ std::string   CAL_DiscrLO;
      _EXTERN_ std::string   CAL_DiscrHI;

      //   ... more for CAL calib types

      //   ... more for ACD calib types 

      //   ... simple type whose "data" just come from MySQL metadata row
      //       for testing
      _EXTERN_ std::string   TestMetadataInfo;

      typedef  std::vector<std::pair <std::string, CLID> > CalibPairCol;
      typedef  CalibPairCol::const_iterator PairIt;
      _EXTERN_    CalibPairCol pairs;
    }
    

#undef _EXTERN_
#endif 
