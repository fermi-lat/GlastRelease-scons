// $Header$

#define _CalibData_CalibModel_cxx

#include "CalibData/CalibModel.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernl/ClassID.h"

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

class calibModel {
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

    CalibData::TKR_MIPEfficiency = tkr + "_MIPEff";

    CalibData::CAL_LightAtt = cal + "LightAtt";
    CalibData::CAL_LightAsym = cal + "LightAsym";
    CalibData::CAL_LightYield = cal + "LightYield";
    CalibData::CAL_Ped = cal + "Ped";
    CalibData::CAL_ElecGain = cal + "ElectGain";
    CalibData::CAL_IntNonlin = cal + "IntNonlin";
    CalibData::CAL_DiffNonlin = cal + "DiffNonlin";
    CalibData::CAL_HotChan = cal + "HotChan";
    CalibData::CAL_DeadChan = cal + "DeadChan";
    CalibData::CAL_DiscrLO = cal + "DiscrLO";
    CalibData::CAL_DiscrHI = cal + "DiscrHI";

    CalibData::TestMetadataInfo = test + "_MetadataInfo";

    CalibData::CalibPairCol.push_back(std::pair(CalibData::TKR_BadChan,
                                                CLID_Calib_TKR_BadChan));
    // Use same class for hot strips, dead strips or (merged) bad strips,
    // but different path in TDDS
    CalibData::CalibPairCol.push_back(std::pair(CalibData::TKR_HotChan,
                                                CLID_Calib_TKR_BadChan));
    CalibData::CalibPairCol.push_back(std::pair(CalibData::TKR_DeadChan,
                                                CLID_Calib_TKR_BadChan));
    CalibData::CalibPairCol.push_back(std::pair(CalibData::TestMetadataInfo,
                                                CLID_Calib_MetadataInfo));
    /* or maybe..
    CalibData::CalibPairCol.push_back(std::make_pair(CalibData::TkrBadStrips,
                                                CLID_Calib_TkrBadStrips));
    CalibData::CalibPairCol.push_back(std::make_pair(CalibData::TestMetadataInfo,
                                                CLID_Calib_MetadataInfo));
    */

  }

};       // end of calibModel class definition

// Instantiate an instance to get the ball rolling.
static CalibModel mod;

// Start class ids at 6000 to stay well away from Gaudi classes and
// our TDS event classes.  Layout could be
//     6000 - 6099        internals
//     6100 - 6199        tracker
//     6200 - 6299        calorimeter
//     6300 - 6399        ACD
//     6400 - 6499        this space intentionally left blank in case
//                        we have calibrations spanning subsystems
//     6500 - 6599        test 

const CLID& CLID_Calib_CalibCLIDNode     = 6000;
const CLID& CLID_Calib_CalibBase         = 6001;
const CLID& CLID_Calib_TkrBadStrips      = 6100;

const CLID& CLID_Calib_MetadataInfo       = 6500;

#undef _CalibData_CalibModel_cxx




