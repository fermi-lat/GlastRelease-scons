// $Header$
#ifndef CalibData_AncCalibTaggerPed_h
#define CalibData_AncCalibTaggerPed_h

#include "CalibData/Anc/AncCalibBase.h"
#include "CalibData/Anc/AncTaggerPed.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  class AncCalibTaggerPed : public AncCalibBase {

  public:
    AncCalibTaggerPed(unsigned nMod=4, unsigned nLay=2, unsigned nChan=384);

    ~AncCalibTaggerPed();

    /// Override putRange implementation in order to add consistency
    /// check
    bool putChan(unsigned iMod, unsigned iLay, unsigned iChan, 
                 RangeBase* data);

    virtual const CLID& clID() const {return classID(); }

    static const CLID& classID();

  };
}
#endif
