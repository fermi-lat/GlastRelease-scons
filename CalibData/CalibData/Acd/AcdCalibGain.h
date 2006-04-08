// $Header$
#ifndef CalibData_AcdCalibGain_h
#define CalibData_AcdCalibGain_h

#include "CalibData/Acd/AcdCalibBase.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  class AcdCalibGain : public AcdCalibBase {

  public:
    AcdCalibGain(unsigned nFace=5, unsigned nRow=5, unsigned nCol=5, 
                 unsigned nNA=11, unsigned nPmt=2);

    ~AcdCalibGain();

    /// Override putRange implementation in order to add consistency
    /// check
    bool putPmt(idents::AcdId id, unsigned pmt, RangeBase* data);

    bool putPmt(unsigned face, unsigned row, unsigned col, 
                unsigned pmt, RangeBase* data);

    virtual const CLID& clID() const {return classID(); }

    static const CLID& classID();

  };
}
#endif
