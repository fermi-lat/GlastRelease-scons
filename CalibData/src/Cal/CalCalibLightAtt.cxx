// $Header$

#include "CalibData/Cal/CalCalibLightAtt.h"
#include "CalibData/Cal/CalFinder.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  CalCalibLightAtt::CalCalibLightAtt(unsigned nTowerRow, unsigned nTowerCol, 
                           unsigned nLayer, unsigned nXtal, 
                           unsigned nFace, unsigned nRange) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, nFace, nRange) {
    unsigned ix = 0;
    unsigned size = m_finder->getSize();

    LightAtt* pLightAtts = new LightAtt[size];
    for (ix = 0; ix < size; ix++) {
      m_ranges[ix] = pLightAtts; 
      ++pLightAtts;
    }
  }

  CalCalibLightAtt::~CalCalibLightAtt() {
    LightAtt* pLightAtts = dynamic_cast<LightAtt*>(m_ranges[0]);
    if (pLightAtts) delete [] pLightAtts;
  }

  const CLID& CalCalibLightAtt::classID()   {return CLID_Calib_CAL_LightAtt;}

  bool CalCalibLightAtt::putRange(idents::CalXtalId id, unsigned range, 
                             unsigned face, RangeBase* data) {
    if (!dynamic_cast<LightAtt* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }

  bool CalCalibLightAtt::putRange(unsigned towerRow, unsigned towerCol, 
                             unsigned layer, unsigned xtal, unsigned range,
                             unsigned face, RangeBase* data) {
    if (!dynamic_cast<LightAtt* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                  range, face, data);
  }

}
