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

std::string root;
std::string tkr;
std::string cal;
std::string acd;
std::string test;
std::string nas;
std::string anc;

class CalibModel {
public:
  /** The constructor sets values into the externally-accessible
      string variables */
  CalibModel () {
    // Initialize a bunch of strings here
    // First few are just for convenience in assembling the    
    // the rest; no need for them to be public.
    std::string root = "/Calib";
    std::string acd = root + "/ACD";
    std::string cal = root + "/CAL";
    std::string tkr = root + "/TKR";
    std::string test = root + "/Test";
    std::string nas = root + "/NAS";
    std::string anc = root + "/ANC";

    // Practically all of these don't exist in any shape or form,
    // but go ahead an reserve names anyway.

    // An alternative way to do this would be to have some appropriate
    // service (CalibDataSvc?  ICalibMetadataCnvSvc?) look up all
    // possible names known to calibUtil.
    CalibData::ACD_Eff = acd +"_Eff";
    CalibData::ACD_ThreshHigh = acd +"_ThreshHigh";
    CalibData::ACD_ThreshVeto = acd +"_ThreshVeto";
    CalibData::ACD_Ped = acd +"_Ped";
    CalibData::ACD_ElecGain = acd +"_ElecGain";
    CalibData::TKR_BadChan = tkr + "_BadChan";
    CalibData::TKR_HotChan = tkr + "_HotChan";
    CalibData::TKR_DeadChan = tkr + "_DeadChan";

    CalibData::TKR_TOTSignal = tkr +"_TOTSignal";
    CalibData::TKR_TOTDist = tkr + "_TOTDist";

    CalibData::TKR_MIPEff = tkr + "_MIPEff";

    CalibData::TKR_Splits = tkr + "_Splits";
  
    CalibData::TKR_ChargeScale = tkr + "_ChargeScale";
    CalibData::TKR_TrgThresh = tkr + "_TrgThresh";
    CalibData::TKR_DataThresh = tkr + "_DataThresh";

    CalibData::CAL_LightAtt = cal + "_LightAtt";
    CalibData::CAL_LightAsym = cal + "_LightAsym";
    CalibData::CAL_LightYield = cal + "_LightYield";
    CalibData::CAL_Ped = cal + "_Ped";
    CalibData::CAL_ElecGain = cal + "_ElecGain";
    CalibData::CAL_IntNonlin = cal + "_IntNonlin";
    CalibData::CAL_DiffNonlin = cal + "_DiffNonlin";
    CalibData::CAL_HotChan = cal + "_HotChan";
    CalibData::CAL_DeadChan = cal + "_DeadChan";

    // As of Sept 2004, CAL folk expect never to use DiscrLO and DiscrHI
    //    CalibData::CAL_DiscrLO = cal + "_DiscrLO";
    //    CalibData::CAL_DiscrHI = cal + "_DiscrHI";
    CalibData::CAL_MuSlope = cal + "_MuSlope";

    CalibData::CAL_MevPerDac = cal + "_MevPerDac";
    CalibData::CAL_TholdCI   = cal + "_TholdCI";
    CalibData::CAL_TholdMuon = cal + "_TholdMuon";
    CalibData::CAL_Asym      = cal + "_Asym";

    CalibData::NAS_TowerCfg = nas + "_TowerCfg";

    CalibData::ANC_TaggerPed = anc + "_TaggerPed";
    CalibData::ANC_TaggerGain = anc + "_TaggerGain";
    CalibData::ANC_QdcPed = anc + "_QdcPed";

    CalibData::Test_Gen = test + "_Gen";
    CalibData::Test_1   = test + "_1";


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


    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_Ped,
                                              CLID_Calib_ACD_Ped));
    CalibData::pairs.push_back(std::make_pair(CalibData::ACD_ElecGain,
                                              CLID_Calib_ACD_ElecGain));

    CalibData::pairs.push_back(std::make_pair(CalibData::Test_1,
                                              CLID_Calib_CalibTest1));

    CalibData::pairs.push_back(std::make_pair(CalibData::NAS_TowerCfg,
                                              CLID_Calib_NAS_TowerCfg));

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




