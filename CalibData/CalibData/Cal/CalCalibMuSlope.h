// $Header$
#ifndef CalibData_CalCalibMuSlope_h
#define CalibData_CalCalibMuSlope_h

#include "CalibData/Cal/CalCalibBase.h"
#include "CalibData/Cal/MuSlope.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  class CalCalibMuSlope : public CalCalibBase {

  public:
    CalCalibMuSlope(unsigned nTowerRow=4, unsigned nTowerCol=4, 
                    unsigned nLayer=8, 
                    unsigned nXtal=12, unsigned nFace=1, unsigned nRange=4);

    ~CalCalibMuSlope();

    /// Override putRange implementation in order to add consistency
    /// check
    bool putRange(idents::CalXtalId id, unsigned range, 
                  unsigned face, RangeBase* data);

    bool putRange(unsigned towerRow, unsigned towerCol, 
                  unsigned layer, unsigned xtal, unsigned range,
                  unsigned face, RangeBase* data);

    virtual const CLID& clID() const {return classID(); }

    static const CLID& classID();


  };

}
#endif
