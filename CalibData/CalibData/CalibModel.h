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


// Start class ids at 6000 to stay well away from Gaudi classes and
// our TDS event classes.  Layout could be
//     6000 - 6099        internals
//     6100 - 6199        tracker
//     6200 - 6299        calorimeter
//     6300 - 6399        ACD
//     6400 - 6499        this space intentionally left blank in case
//                        we have calibrations spanning subsystems
//     6500 - 6599        test 
//     6600 - 6699        NAS: anything not associated with a single subsystem
//                        

const CLID CLID_Calib_CalibCLIDNode     = 6000;
const CLID CLID_Calib_CalibBase         = 6001;

const CLID CLID_Calib_TKR_HotChan       = 6100;
const CLID CLID_Calib_TKR_DeadChan      = 6101;
const CLID CLID_Calib_TKR_BadChan      = 6102;
const CLID CLID_Calib_TKR_TOTSignal      = 6103;
const CLID CLID_Calib_TKR_TOTDist      = 6104;
const CLID CLID_Calib_TKR_MIPEff       = 6105;
const CLID CLID_Calib_TKR_Splits       = 6106;
const CLID CLID_Calib_TKR_ChargeScale  = 6107;
const CLID CLID_Calib_TKR_TrgThresh    = 6108;
const CLID CLID_Calib_TKR_DataThresh   = 6109;


const CLID CLID_Calib_CAL_LightAtt   = 6200;
const CLID CLID_Calib_CAL_LightAsym  = 6201;
const CLID CLID_Calib_CAL_LightYield = 6202;
const CLID CLID_Calib_CAL_Ped        = 6203;
const CLID CLID_Calib_CAL_ElecGain   = 6204;
const CLID CLID_Calib_CAL_IntNonlin  = 6205;
const CLID CLID_Calib_CAL_DiffNonlin = 6206;
const CLID CLID_Calib_CAL_HotChan    = 6207;
const CLID CLID_Calib_CAL_DeadChan   = 6208;
// As of Sept 2004, CAL folk expect never to use DiscrLO and DiscrHI
// const CLID CLID_Calib_CAL_DiscrLO    = 6209;
// const CLID CLID_Calib_CAL_DiscrHI    = 6210;
const CLID CLID_Calib_CAL_MuSlope    = 6211;

//  Added 13 Sept. 2004
const CLID CLID_Calib_CAL_MevPerDac  = 6212;
const CLID CLID_Calib_CAL_TholdCI    = 6213;
const CLID CLID_Calib_CAL_TholdMuon  = 6214;
const CLID CLID_Calib_CAL_Asym       = 6215;

const CLID CLID_Calib_ACD_Eff        = 6300;
const CLID CLID_Calib_ACD_ThreshHigh = 6301;
const CLID CLID_Calib_ACD_ThreshVeto = 6302;
const CLID CLID_Calib_ACD_Ped        = 6304;
const CLID CLID_Calib_ACD_ElecGain   = 6305;

const CLID CLID_Calib_CalibTest1       = 6500;

const CLID CLID_Calib_NAS_TowerCfg     = 6600;



// For everybody except the CalibModel class implementation file,
// the variables are extern.  CalibModel.cxx actually defines them.

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
      _EXTERN_ std::string   TKR_Splits;

      _EXTERN_ std::string   TKR_ChargeScale;
      _EXTERN_ std::string   TKR_TrgThresh;
      _EXTERN_ std::string   TKR_DataThresh;
    

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
      _EXTERN_ std::string   CAL_MuSlope;

      _EXTERN_ std::string   CAL_MevPerDac;
      _EXTERN_ std::string   CAL_TholdCI;
      _EXTERN_ std::string   CAL_TholdMuon;
      _EXTERN_ std::string   CAL_Asym;
      //   ... more for CAL calib types

      //   ... more for ACD calib types 

      //   ... simple type whose "data" just come from MySQL metadata row
      //       for testing
      _EXTERN_ std::string   Test_Gen;

      // Simple xml test type
      _EXTERN_ std::string   Test_1;

      // cross-subsystem types
      _EXTERN_ std::string   NAS_TowerCfg;
      typedef  std::vector<std::pair <std::string, CLID> > CalibPairCol;
      typedef  CalibPairCol::const_iterator PairIt;
      _EXTERN_    CalibPairCol pairs;
    }
    

#undef _EXTERN_
#endif 
