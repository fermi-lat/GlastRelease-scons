// $Header$
#ifndef CalibData_CalCalibLightAsym_h
#define CalibData_CalCalibLightAsym_h

#include "CalibData/Cal/CalCalibBase.h"
#include "CalibData/Cal/LightAsym.h"
#include "CalibData/CalibModel.h"

namespace CalibData {
  class RangeBase;

  class CalCalibLightAsym : public CalCalibBase {

  public:
    CalCalibLightAsym(unsigned nTowerRow=4, unsigned nTowerCol=4, 
                 unsigned nLayer=8, 
                 unsigned nXtal=12, unsigned nFace=1, unsigned nRange=2);

    ~CalCalibLightAsym();

    /// Override so that both LEX1 and LEX8 range are associated with
    /// same data and similarly for HEX1, HEX8
    RangeBase* getRange(idents::CalXtalId id, unsigned range=0,
                        unsigned face = 0);

    /// Override putRange implementation in order to add consistency
    /// check; also to translate "range" correctly 
    bool putRange(idents::CalXtalId id, unsigned range, 
                  unsigned face, RangeBase* data);

    bool putRange(unsigned towerRow, unsigned towerCol, 
                  unsigned layer, unsigned xtal, unsigned range,
                  unsigned face, RangeBase* data);

    virtual const CLID& clID() const {return classID(); }

    static const CLID& classID();

  private:
    unsigned m_nRange;
  };

}
#endif
