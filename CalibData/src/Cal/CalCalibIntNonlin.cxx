// $Header$

#include "CalibData/Cal/CalCalibIntNonlin.h"
#include "CalibData/Cal/CalFinder.h"
#include "CalibData/DacCol.h"
#include "CalibData/CalibModel.h"

namespace CalibData {

  CalCalibIntNonlin::CalCalibIntNonlin(unsigned nTowerRow, unsigned nTowerCol, 
                                       unsigned nLayer, unsigned nXtal, 
                                       unsigned nFace, unsigned nRange,
                                       unsigned nDacCol) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, nFace, nRange, nDacCol)
  {
      cGuts(nTowerRow, nTowerCol, nLayer, nXtal,
            nFace, nRange, nDacCol);
  }

  CalCalibIntNonlin::~CalCalibIntNonlin() {
    IntNonlin* pIntNonlins = dynamic_cast<IntNonlin* >(m_ranges[0]);
    if (pIntNonlins) delete [] pIntNonlins;
    if (m_dacCols.size() > 0) {
      DacCol* pDacCols = dynamic_cast<DacCol* > (m_dacCols[0]);
      if (pDacCols) delete [] pDacCols;  
    }
  }

  const CLID& CalCalibIntNonlin::classID()   {return CLID_Calib_CAL_IntNonlin;}


  bool CalCalibIntNonlin::putRange(idents::CalXtalId id, unsigned range, 
                             unsigned face, RangeBase* data) {
    if (!dynamic_cast<IntNonlin* >(data)) return false;
    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }

  bool CalCalibIntNonlin::putRange(unsigned towerRow, unsigned towerCol, 
                              unsigned layer, unsigned xtal, unsigned range,
                              unsigned face, RangeBase* data) {
    if (!dynamic_cast<IntNonlin* >(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                  range, face, data);

  }

  void CalCalibIntNonlin::cGuts(unsigned /* nTowerRow */, 
                                unsigned /* nTowerCol */, 
                                unsigned /* nLayer */, unsigned /* nXtal */, 
                                unsigned /* nFace */, unsigned /* nRange */, 
                                unsigned nDacCol) {

    unsigned ix = 0;
    unsigned size = m_finder->getSize();
    
    IntNonlin* pIntNonlins = new IntNonlin[size];
    for (ix = 0; ix < size; ix++) {
      m_ranges[ix] = pIntNonlins; 
      ++pIntNonlins;
    }
    if (nDacCol != 0) {
      DacCol* pDacCols = new DacCol[nDacCol];
      for (ix = 0; ix < nDacCol; ix++) {
        m_dacCols.push_back(pDacCols);
        ++pDacCols;
      }
    }
  }

}
