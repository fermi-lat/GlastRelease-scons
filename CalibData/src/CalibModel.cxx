// $Header$

#define _CalibData_CalibModel_cxx

#include "CalibData/CalibModel.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ClassID.h"

/**  @file CalibModel.cxx 
     Implementation for CalibModel class, which initializes strings
     for use as paths for TDDS data.  Also includes definitions
     of Gaudi class ids for TDDS DataObjects.
*/

/**  @class CalibModel
     Provides convenience definitions of strings which are paths
     to calibration data in the TDDS.  The class has no data members.
     The strings it initializes are just "there"; existing in 
     this module and extern everywhere else (see the fancy footwork
     in CalibModel.h).
     The single static instance of CalibModel declared in this file
     causes the constructor to be run, initializing the strings, and
     a vector of pairs.  The pairs are used to match each string
     with the appropriate Gaudi class id.

     NOTE:  This scheme might have to be revised.  We haven't dealt
            with different instruments, nor with flavors.  Assuming,
            as seems likely, that we *never* have to keep track of
            constants for more than one instrument simultaneously 
            (most likely not even sequentially within a single job),
            we can leave any mention of instrument out of these
            paths.  If we want to preserve the possibility of maintaining
            more than one flavor simultaneously, flavor will have to
            be part of the path string.  We could tack it on the end
            dynamically without having to do any violence to what
            we've got below.
*/

class CalibModel {
public:
  /** The constructor sets values into the externally-accessible
      string variables */
  CalibModel () {
    // Initialize a bunch of strings here
    // First few are just for convenience in assembling the    
    // the rest
    CalibData::root = "/Calib";
    CalibData::acd = CalibData::root + "/ACD";
    CalibData::cal = CalibData::root + "/CAL";
    CalibData::tkr = CalibData::root + "/TKR";
    CalibData::test = CalibData::root + "/Test";
    CalibData::nas = CalibData::root + "/NAS";
    CalibData::anc = CalibData::root + "/ANC";

    // Practically all of these don't exist in any shape or form,
    // but go ahead an reserve names anyway.

    // An alternative way to do this would be to have some appropriate
    // service (CalibDataSvc?  ICalibMetadataCnvSvc?) look up all
    // possible names known to calibUtil.
    CalibData::ACD_Eff = CalibData::acd +"_Eff";
    CalibData::ACD_ThreshHigh = CalibData::acd +"_ThreshHigh";
    CalibData::ACD_ThreshVeto = CalibData::acd +"_ThreshVeto";
    CalibData::ACD_Ped = CalibData::acd +"_Ped";
    CalibData::ACD_ElecGain = CalibData::acd +"_ElecGain";
    CalibData::ACD_Range = CalibData::acd +"_Range";
    CalibData::ACD_HighRange = CalibData::acd +"_HighRange";
    CalibData::ACD_CoherentNoise = CalibData::acd +"_CoherentNoise";
    CalibData::ACD_Ribbon = CalibData::acd +"_Ribbon";
    CalibData::ACD_HighPed = CalibData::acd +"_HighPed";
    CalibData::ACD_Carbon = CalibData::acd +"_Carbon";
    CalibData::ACD_VetoFit = CalibData::acd +"_VetoFit";
    CalibData::ACD_CnoFit = CalibData::acd +"_CnoFit";

    CalibData::TKR_BadChan = CalibData::tkr + "_BadChan";
    CalibData::TKR_HotChan = CalibData::tkr + "_HotChan";
    CalibData::TKR_DeadChan = CalibData::tkr + "_DeadChan";

    CalibData::TKR_TOTSignal = CalibData::tkr +"_TOTSignal";
    CalibData::TKR_TOTDist = CalibData::tkr + "_TOTDist";

    CalibData::TKR_MIPEff = CalibData::tkr + "_MIPEff";

    CalibData::TKR_Splits = CalibData::tkr + "_Splits";
  
    CalibData::TKR_ChargeScale = CalibData::tkr + "_ChargeScale";
    CalibData::TKR_TrgThresh = CalibData::tkr + "_TrgThresh";
    CalibData::TKR_DataThresh = CalibData::tkr + "_DataThresh";
    CalibData::TKR_TowerAlign = CalibData::tkr + "_TowerAlign";
    CalibData::TKR_InternalAlign = CalibData::tkr + "_InternalAlign";

    CalibData::CAL_LightAtt = CalibData::cal + "_LightAtt";
    CalibData::CAL_LightAsym = CalibData::cal + "_LightAsym";
    CalibData::CAL_LightYield = CalibData::cal + "_LightYield";
    CalibData::CAL_Ped = CalibData::cal + "_Ped";
    CalibData::CAL_ElecGain = CalibData::cal + "_ElecGain";
    CalibData::CAL_IntNonlin = CalibData::cal + "_IntNonlin";
    CalibData::CAL_DiffNonlin = CalibData::cal + "_DiffNonlin";
    CalibData::CAL_HotChan = CalibData::cal + "_HotChan";
    CalibData::CAL_DeadChan = CalibData::cal + "_DeadChan";

    // As of Sept 2004, CAL folk expect never to use DiscrLO and DiscrHI
    //    CalibData::CAL_DiscrLO = cal + "_DiscrLO";
    //    CalibData::CAL_DiscrHI = cal + "_DiscrHI";
    CalibData::CAL_MuSlope = CalibData::cal + "_MuSlope";

    CalibData::CAL_MevPerDac = CalibData::cal + "_MevPerDac";
    CalibData::CAL_TholdCI   = CalibData::cal + "_TholdCI";
    CalibData::CAL_TholdMuon = CalibData::cal + "_TholdMuon";
    CalibData::CAL_Asym      = CalibData::cal + "_Asym";

    CalibData::NAS_TowerCfg    = CalibData::nas + "_TowerCfg";
    CalibData::NAS_SAABoundary = CalibData::nas + "_SAABoundary";
    CalibData::NAS_LATAlignment = CalibData::nas + "_LATAlignment";

    CalibData::ANC_TaggerPed = CalibData::anc + "_TaggerPed";
    CalibData::ANC_TaggerGain = CalibData::anc + "_TaggerGain";
    CalibData::ANC_QdcPed = CalibData::anc + "_QdcPed";

    CalibData::Test_Gen = CalibData::test + "_Gen";
    CalibData::Test_1   = CalibData::test + "_1";


    // Use same class for hot strips, dead strips or (merged) bad strips,
    // but different path in TDDS

    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_HotChan,
                                              CLID_Calib_TKR_BadChan));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_DeadChan,
                                              CLID_Calib_TKR_BadChan));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_Splits,
                                              CLID_Calib_TKR_Splits));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_TOTSignal,
                                              CLID_Calib_TKR_TOTSignal));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_ChargeScale,
                                              CLID_Calib_TKR_ChargeScale));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_TrgThresh,
                                              CLID_Calib_TKR_TrgThresh));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_DataThresh,
                                              CLID_Calib_TKR_DataThresh));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_TowerAlign,
                                              CLID_Calib_TKR_TowerAlign));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_InternalAlign,
                                              CLID_Calib_TKR_InternalAlign));




    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_Ped,
                                              CLID_Calib_CAL_Ped));
    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_ElecGain,
                                              CLID_Calib_CAL_ElecGain));
    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_MuSlope,
                                              CLID_Calib_CAL_MuSlope));
    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_LightAtt,
                                              CLID_Calib_CAL_LightAtt));
    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_LightAsym,
                                              CLID_Calib_CAL_LightAsym));
    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_IntNonlin,
                                              CLID_Calib_CAL_IntNonlin));

    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_MevPerDac,
                                              CLID_Calib_CAL_MevPerDac));
    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_TholdCI,
                                              CLID_Calib_CAL_TholdCI));
    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_TholdMuon,
                                              CLID_Calib_CAL_TholdMuon));
    CalibData::pairs.push_back(std::make_pair(CalibData::CAL_Asym,
                                              CLID_Calib_CAL_Asym));


    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_ThreshVeto,
                                              CLID_Calib_ACD_ThreshVeto));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_ThreshHigh,
                                              CLID_Calib_ACD_ThreshHigh));

    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_ElecGain,
                                              CLID_Calib_ACD_ElecGain));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_Ped,
                                              CLID_Calib_ACD_Ped));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_Range,
                                              CLID_Calib_ACD_Range));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_HighRange,
                                              CLID_Calib_ACD_HighRange));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_CoherentNoise,
                                              CLID_Calib_ACD_CoherentNoise));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_Ribbon,
                                              CLID_Calib_ACD_Ribbon));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_HighPed,
                                              CLID_Calib_ACD_HighPed));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_Carbon,
                                              CLID_Calib_ACD_Carbon));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_VetoFit,
                                              CLID_Calib_ACD_VetoFit));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_CnoFit,
                                              CLID_Calib_ACD_CnoFit));

    CalibData::pairs.push_back(std::make_pair(CalibData::Test_1,
                                              CLID_Calib_CalibTest1));

    CalibData::pairs.push_back(std::make_pair(CalibData::NAS_TowerCfg,
                                              CLID_Calib_NAS_TowerCfg));
    CalibData::pairs.push_back(std::make_pair(CalibData::NAS_SAABoundary,
                                              CLID_Calib_NAS_SAABoundary));
    CalibData::pairs.push_back(std::make_pair(CalibData::NAS_LATAlignment,
                                              CLID_Calib_NAS_LATAlignment));

    CalibData::pairs.push_back(std::make_pair(CalibData::ANC_TaggerPed,
                                              CLID_Calib_ANC_TaggerPed));

    CalibData::pairs.push_back(std::make_pair(CalibData::ANC_TaggerGain,
                                              CLID_Calib_ANC_TaggerGain));

    CalibData::pairs.push_back(std::make_pair(CalibData::ANC_QdcPed,
                                              CLID_Calib_ANC_QdcPed));


  }

};       // end of calibModel class definition

// Instantiate an instance to get the ball rolling.
static CalibModel mod;


#undef _CalibData_CalibModel_cxx




