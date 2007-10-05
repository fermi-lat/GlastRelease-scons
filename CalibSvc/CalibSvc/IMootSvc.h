//$Header$
#ifndef IMootSvc_h
#define IMootSvc_h 1

#include "GaudiKernel/IInterface.h"

// External constants
#include "GaudiKernel/ClassID.h"


// Forward declarations
class DataObject;

namespace MOOT {
  class MootQuery;
}

static const InterfaceID IID_IMootSvc ("IMootSvc", 1, 0);

/** @class IMootSvc 

    Abstract interface of a conversion service for GLAST config meta data
    persistency (to start; may add calibration later)

    Handles creation and updating condition data objects (i.e. DataObjects
    implementing IValidity).

    Adapted from Andrea Valassi's IConditionsDBCnvSvc interface

    @author Joanne Bogart
    @date September 2007
*/

class IMootSvc : virtual public IInterface   {

public:
  // Re-implemented from IInterface
  static const InterfaceID& interfaceID() { return IID_IMootSvc; }

  // Start with non-conversion service stuff
  // Get TDS path for LATC source files
  virtual std::string getLatcSourcePath()=0;

  /// Return index in MootParmCol of specified class
  virtual int latcParmIx(const std::string& parmClass)=0;

  // Enumeration of all precincts which may have associated parameter classes.
  // Does not include "container precints" like LPA or LCI_GLOBAL_ACD
  typedef enum {
    PR_generic = 0,
    PR_ACD_Mode,
    PR_firstLatc = PR_ACD_Mode,
    PR_ACD_Bias,
    PR_ACD_Hld,
    PR_ACD_PHA,
    PR_ACD_Veto,
    PR_ACD_Timing,
    PR_CAL_Mode,
    PR_CAL_Timing,
    PR_CAL_LAC,
    PR_CAL_FLE,
    PR_CAL_FHE,
    PR_CAL_ULD,
    PR_TKR_Mode,
    PR_TKR_Timing,
    PR_TKR_Strips,
    PR_TKR_Thresh,
    PR_GNL_Mode,
    PR_GNL_Timing,
    PR_TRG_ROI,
    PR_TRG_GEM,
    PR_lastLatc = PR_TRG_GEM,
    PR_ACD_LCI,
    PR_CAL_LCI,
    PR_TKR_LCI,
    PR_count
  }  Precincts;

  // Given a (latc??) precinct, return  list of indices? iterators?



 public:

  

  

  // Get handle for metadata access from mootCore.
  virtual MOOT::MootQuery* getConnection(bool verbose=true) = 0;


};

#endif

