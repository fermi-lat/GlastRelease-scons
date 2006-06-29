// $Header$
#ifndef CalibData_AncCalibQdcPed_h
#define CalibData_AncCalibQdcPed_h

#include "CalibData/Anc/AncCalibBase.h"
#include "CalibData/Anc/AncQdcPed.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  class AncCalibQdcPed : public AncCalibBase {

  public:
    AncCalibQdcPed(unsigned nMod=2, unsigned nChan=32);

    ~AncCalibQdcPed();

    /// Override putRange implementation in order to add consistency
    /// check
    bool putChan(unsigned iMod, unsigned iLay, unsigned iChan, 
                 RangeBase* data);

    /// Since qdc channels are identified by module and channel only
    /// (no layer) define this additional 'get' method so apps have
    /// a natural way to access channel data.
    RangeBase* getQdcChan(unsigned iMod=0, unsigned iChan=0) {
      return getChan(iMod, 0, iChan);
    }

    virtual const CLID& clID() const {return classID(); }

    static const CLID& classID();

  };
}
#endif
