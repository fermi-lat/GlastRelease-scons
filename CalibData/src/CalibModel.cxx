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

    CalibData::CAL_LightAtt = cal + "_LightAtt";
    CalibData::CAL_LightAsym = cal + "_LightAsym";
    CalibData::CAL_LightYield = cal + "_LightYield";
    CalibData::CAL_Ped = cal + "_Ped";
    CalibData::CAL_ElecGain = cal + "_ElectGain";
    CalibData::CAL_IntNonlin = cal + "_IntNonlin";
    CalibData::CAL_DiffNonlin = cal + "_DiffNonlin";
    CalibData::CAL_HotChan = cal + "_HotChan";
    CalibData::CAL_DeadChan = cal + "_DeadChan";
    CalibData::CAL_DiscrLO = cal + "_DiscrLO";
    CalibData::CAL_DiscrHI = cal + "_DiscrHI";

    CalibData::Test_Gen = test + "_Gen";


    // Use same class for hot strips, dead strips or (merged) bad strips,
    // but different path in TDDS

    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_HotChan,
                                                CLID_Calib_TKR_HotChan));
    CalibData::pairs.push_back(std::make_pair(CalibData::TKR_DeadChan,
                                                CLID_Calib_TKR_DeadChan));


  }

};       // end of calibModel class definition

/*
const CLID& CLID_Calib_CalibCLIDNode     = 6000;
const CLID& CLID_Calib_CalibBase         = 6001;

const CLID& CLID_Calib_TKR_HotChan       = 6100;
const CLID& CLID_Calib_TKR_DeadChan      = 6101;
const CLID& CLID_Calib_TKR_BadChan      = 6102;
const CLID& CLID_Calib_TKR_TOTSignal      = 6103;
const CLID& CLID_Calib_TKR_TOTDist      = 6104;
const CLID& CLID_Calib_TKR_MIPEff       = 6105;

const CLID& CLID_Calib_CAL_LightAtt   = 6200;
const CLID& CLID_Calib_CAL_LightAsym  = 6201;
const CLID& CLID_Calib_CAL_LightYield = 6202;
const CLID& CLID_Calib_CAL_Ped        = 6203;
const CLID& CLID_Calib_CAL_ElecGain   = 6204;
const CLID& CLID_Calib_CAL_IntNonlin  = 6205;
const CLID& CLID_Calib_CAL_DiffNonlin = 6206;
const CLID& CLID_Calib_CAL_HotChan    = 6207;
const CLID& CLID_Calib_CAL_DeadChan   = 6208;
const CLID& CLID_Calib_CAL_DiscrLO    = 6209;
const CLID& CLID_Calib_CAL_DiscrHI    = 6210;

const CLID& CLID_Calib_ACD_Eff        = 6300;
const CLID& CLID_Calib_ACD_ThreshHigh = 6301;
const CLID& CLID_Calib_ACD_ThreshVeto = 6302;
const CLID& CLID_Calib_ACD_Ped        = 6304;
const CLID& CLID_Calib_ACD_ElecGain   = 6305;

const CLID& CLID_Calib_CalibTest1       = 6500;
*/
// Instantiate an instance to get the ball rolling.
static CalibModel mod;


#undef _CalibData_CalibModel_cxx




