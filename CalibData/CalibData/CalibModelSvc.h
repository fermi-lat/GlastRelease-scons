// $Header$
#ifndef CalibData_CalibModelSvc_h
#define CalibData_CalibModelSvc_h

/**
   @file CalibModelSvc.h

  A nearly-contentless class which only exists to make information
  from CalibModel available to components other than CalibData.

  For now just serve information needed by CalibDataSvc.
*/
#include <string>
#include <vector>
#include "GaudiKernel/ClassID.h"

namespace CalibData {
  class CalibModelSvc {
  public:
    typedef std::pair<std::string, CLID> CalibPair;

    const std::vector<CalibPair>& getPairs() const;

    CLID getCLIDNodeCLID() const;

    /// Return flavor string, extracted from full path name in TCDS
    static std::string getFlavor(const std::string& fullpath);

    /// Return calibration type name, extracted from full path name in TCDS
    static std::string getCalibType(const std::string& fullpath);
  };
}







#endif
