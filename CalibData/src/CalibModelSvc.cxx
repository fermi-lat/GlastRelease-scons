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


  std::string CalibModelSvc::getFlavor(const std::string& fullpath) {
    // flavor is the last field
    unsigned lastSlash = fullpath.rfind("/");
    return fullpath.substr(lastSlash+1, fullpath.size() - lastSlash);
  }

  std::string CalibModelSvc::getCalibType(const std::string& fullpath) {
    // flavor is the last field
    unsigned lastSlash = fullpath.rfind("/");

    // calib type is second-to-last
    unsigned prevSlash = fullpath.rfind("/", lastSlash - 1);
    return fullpath.substr(prevSlash+1, lastSlash - prevSlash - 1);
  }


}
