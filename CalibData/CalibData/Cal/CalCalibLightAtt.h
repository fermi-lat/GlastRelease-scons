// $Header$
#ifndef CalibData_CalCalibLightAtt_h
#define CalibData_CalCalibLightAtt_h

#include "CalibData/Cal/CalCalibBase.h"
#include "CalibData/Cal/LightAtt.h"


namespace CalibData {

  /**
     @class CalCalibLightAtt

     This class is a container for per-range light attenuation data.  It
     inherits all its implementation from CalCalibBase except for
     knowledge of its own CLID and a check in putRange to make
     sure the data being added is of the right sort.
  */
  class CalCalibLightAtt : public CalCalibBase {

  public:
    CalCalibLightAtt(unsigned nTowerRow=4, unsigned nTowerCol=4, unsigned nLayer=8, 
                unsigned nXtal=12, unsigned nFace=1, unsigned nRange=1);

    ~CalCalibLightAtt();

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
