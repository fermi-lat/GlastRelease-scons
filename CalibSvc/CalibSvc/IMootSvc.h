//$Header$
#ifndef IMootSvc_h
#define IMootSvc_h 1

#include "GaudiKernel/IInterface.h"

// External constants
#include "GaudiKernel/ClassID.h"


// Forward declarations

namespace MOOT {
  class MootQuery;
}

namespace CalibData {
  class MootParm;
  class MootParmCol;
}

static const InterfaceID IID_IMootSvc ("IMootSvc", 1, 0);

/** @class IMootSvc 

    Abstract interface of a service for access to MOOT information.
    See also data class definitions in CalibData/Moot


    @author Joanne Bogart
*/

class IMootSvc : virtual public IInterface   {

public:
  // Re-implemented from IInterface
  static const InterfaceID& interfaceID() { return IID_IMootSvc; }

  // Return pointer to Moot parameter collection.  Also set output
  // arg. hw to current hw key
  virtual const CalibData::MootParmCol* getMootParmCol(unsigned& hw)=0;

  /// Return last LATC master key seen in data
  virtual unsigned getHardwareKey()=0;


  /// Return index in MootParmCol of specified class
  virtual int latcParmIx(const std::string& parmClass) const =0;


  // Given a (latc??) precinct, return  list of indices? iterators?



  // Get handle for metadata access from mootCore.
  virtual MOOT::MootQuery* getConnection() const = 0;



};

#endif

