// $Header$

#include "CalibData/Cal/CalCalibGain.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  const CLID& CalCalibGain::classID()   {return CLID_Calib_CAL_ElecGain;}

  bool CalCalibGain::putRange(idents::CalXtalId id, unsigned range, 
                             unsigned face, RangeBase* data) {
    if (!dynamic_cast<Gain* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }

}
