// $Header$
#ifndef CalibData_CalCalibPed_h
#define CalibData_CalCalibPed_h

#include "CalibData/Cal/CalCalibBase.h"
#include "CalibData/Cal/Ped.h"


namespace CalibData {

  /**
     @class CalCalibPed

     This class is a container for per-range pedestal data.  It
     inherits all its implementation from CalCalibBase except for
     knowledge of its own CLID and a check in putRange to make
     sure the data being added is of the right sort.
  */
  class CalCalibPed : public CalCalibBase {

  public:
    CalCalibPed(unsigned nTowerRow=4, unsigned nTowerCol=4, unsigned nLayer=8, 
                unsigned nXtal=12, unsigned nFace=2, unsigned nRange=4);

    ~CalCalibPed() {}

    /// Override putRange implementations in order to add consistency
    /// check
    bool putRange(idents::CalXtalId id, unsigned range, 
                  unsigned face, RangeBase* data);

    bool putRange(unsigned towerRow, unsigned towerCol, 
                  unsigned layer, unsigned xtal, unsigned range,
                  unsigned face, RangeBase* data);

    virtual const CLID& clID() const { return classID(); }

    static const CLID& classID();
  };

}
#endif
