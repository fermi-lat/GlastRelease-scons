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
    std::string acd = root + "/Acd";
    std::string cal = root + "/Cal";
    std::string tkr = root + "/Tkr";
    std::string test = root + "/Test";

    CalibData::TkrBadStrips = tkr + ".BadStrips";
    CalibData::TestMetadataInfo = test + ".MetadataInfo";

    CalibData::CalibPairCol.push_back(std::pair(CalibData::TkrBadStrips,
                                                CLID_Calib_TkrBadStrips));
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
const CLID& CLID_Calib_TkrBadStrips      = 6100;

const CLID& CLID_Calib_MetadataInfo       = 6500;

#undef _CalibData_CalibModel_cxx




