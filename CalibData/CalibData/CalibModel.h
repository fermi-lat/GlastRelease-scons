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
//     6700 - 6799        ANC: ancillary - not part of the LAT at all
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
const CLID CLID_Calib_TKR_TowerAlign   = 6110;
const CLID CLID_Calib_TKR_InternalAlign = 6111;


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
// 
const CLID CLID_Calib_ACD_Range          = 6306;
const CLID CLID_Calib_ACD_HighRange      = 6307;
const CLID CLID_Calib_ACD_CoherentNoise  = 6308;
const CLID CLID_Calib_ACD_Ribbon         = 6309;

const CLID CLID_Calib_ACD_HighPed        = 6310;
const CLID CLID_Calib_ACD_Carbon         = 6311;
const CLID CLID_Calib_ACD_VetoFit        = 6312;
const CLID CLID_Calib_ACD_CnoFit         = 6313; 
const CLID CLID_Calib_ACD_Merit          = 6314;
const CLID CLID_Calib_ACD_PE             = 6315;

const CLID CLID_Calib_CalibTest1       = 6500;

const CLID CLID_Calib_NAS_TowerCfg     = 6600;
const CLID CLID_Calib_NAS_SAABoundary  = 6610;
const CLID CLID_Calib_NAS_LATAlignment = 6611;

const CLID CLID_Calib_ANC_TaggerPed    = 6700;
const CLID CLID_Calib_ANC_TaggerGain   = 6701;
const CLID CLID_Calib_ANC_QdcPed       = 6702;


// For everybody except the CalibModel class implementation file,
// the variables are extern.  CalibModel.cxx actually defines them.

#if defined(_CalibData_CalibModel_cxx)
#define  _EXTERN_ 
#else
#define  _EXTERN_ extern
#endif

    namespace CalibData {
      // Just for building remainder
      _EXTERN_ std::string root;
      _EXTERN_ std::string acd;
      _EXTERN_ std::string cal;
      _EXTERN_ std::string tkr;
      _EXTERN_ std::string test;
      _EXTERN_ std::string nas;
      _EXTERN_ std::string anc;

      // ACD calib types
      _EXTERN_ std::string   ACD_Eff;
      _EXTERN_ std::string   ACD_ThreshHigh;
      _EXTERN_ std::string   ACD_ThreshVeto;
      _EXTERN_ std::string   ACD_Ped;
      _EXTERN_ std::string   ACD_ElecGain;
      _EXTERN_ std::string   ACD_Range;
      _EXTERN_ std::string   ACD_HighRange;
      _EXTERN_ std::string   ACD_CoherentNoise;
      _EXTERN_ std::string   ACD_Ribbon;
      _EXTERN_ std::string   ACD_HighPed;
      _EXTERN_ std::string   ACD_Carbon;
      _EXTERN_ std::string   ACD_VetoFit;
      _EXTERN_ std::string   ACD_CnoFit;
      _EXTERN_ std::string   ACD_PE;

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


      _EXTERN_ std::string   TKR_TowerAlign;
      _EXTERN_ std::string   TKR_InternalAlign;
    

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

      _EXTERN_ std::string   ANC_TaggerPed;
      _EXTERN_ std::string   ANC_TaggerGain;
      _EXTERN_ std::string   ANC_QdcPed;
      //   ... simple type whose "data" just come from MySQL metadata row
      //       for testing
      _EXTERN_ std::string   Test_Gen;

      // Simple xml test type
      _EXTERN_ std::string   Test_1;

      // cross-subsystem types
      _EXTERN_ std::string   NAS_TowerCfg;
      _EXTERN_ std::string   NAS_SAABoundary;
      _EXTERN_ std::string   NAS_LATAlignment;


      typedef  std::vector<std::pair <std::string, CLID> > CalibPairCol;
      typedef  CalibPairCol::const_iterator PairIt;
      _EXTERN_    CalibPairCol pairs;
    }
    

#undef _EXTERN_
#endif 
