// $Header$
#ifndef CalibData_AcdCalibVeto_h
#define CalibData_AcdCalibVeto_h

#include "CalibData/Acd/AcdCalibBase.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  class AcdCalibVeto : public AcdCalibBase {

  public:
    AcdCalibVeto(unsigned nFace=7, unsigned nRow=5, unsigned nCol=5, 
                 unsigned nNA=11, unsigned nPmt=2);

    ~AcdCalibVeto();

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
