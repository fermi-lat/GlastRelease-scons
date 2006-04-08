// $Header$
#ifndef CalibData_AcdCalibCno_h
#define CalibData_AcdCalibCno_h

#include "CalibData/Acd/AcdCalibBase.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  class AcdCalibCno : public AcdCalibBase {

  public:
    AcdCalibCno(unsigned nFace=7, unsigned nRow=5, unsigned nCol=5, 
                unsigned nNA=11, unsigned nPmt=2);

    ~AcdCalibCno();

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
