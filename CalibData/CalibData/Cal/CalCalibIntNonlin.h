// $Header$
#ifndef CalibData_CalCalibIntNonlin_h
#define CalibData_CalCalibIntNonlin_h

#include "CalibData/Cal/CalCalibBase.h"
#include "CalibData/Cal/IntNonlin.h"
#include "CalibData/CalibModel.h"

namespace CalibData {
  class RangeBase;

  class CalCalibIntNonlin : public CalCalibBase {

  public:
    CalCalibIntNonlin(unsigned nTowerRow=4, unsigned nTowerCol=4, 
                      unsigned nLayer=8, 
                      unsigned nXtal=12, unsigned nFace=2, unsigned nRange=4,
                      unsigned nDacCol=4);

    ~CalCalibIntNonlin();


    /// Override putRange implementation in order to add consistency
    /// check;
    bool putRange(idents::CalXtalId id, unsigned range, 
                  unsigned face, RangeBase* data);

    bool putRange(unsigned towerRow, unsigned towerCol, 
                  unsigned layer, unsigned xtal, unsigned range,
                  unsigned face, RangeBase* data);

    virtual const CLID& clID() const {return classID(); }

    static const CLID& classID();

  private:
    /** Due to bug in gcc, gdb can't find symbols in constructors.  This
        method is called by the constructor and does most of the work
    */
    void cGuts(unsigned nTowerRow, unsigned nTowerCol, 
               unsigned nLayer, unsigned nXtal, unsigned nFace, 
               unsigned nRange, unsigned nDacCol);


  };

}
#endif
