// $Header$

#include "CalibData/CalibModelSvc.h"
#include "CalibData/CalibModel.h"

// extern const CLID& CLID_Calib_CalibCLIDNode;

namespace CalibData {
  const std::vector<CalibModelSvc::CalibPair>& CalibModelSvc::getPairs() 
    const {
    return pairs;
  }

  CLID CalibModelSvc::getCLIDNodeCLID() const {
    return CLID_Calib_CalibCLIDNode;
  }
}
