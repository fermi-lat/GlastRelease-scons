// $Header$

/**
  @file  CalMevPerDac.cxx

  Implementation for new calorimeter gains calibration, consisting of 
  per-crystal value for each size diode, plus per-crystal-end constants
  for Big/Small ratio

  Implementation for two classes may be found here.  @a CalMevPerDac is a class
  for per-crystal information.  @a CalMevPerDacCol contains CalMevPerDac 
  information for the entire detector.
*/
#include "CalibData/Cal/CalMevPerDac.h"
#include "CalibData/Cal/Xpos.h"

namespace CalibData {
  // CalMevPerDac implementation
  CalMevPerDac::CalMevPerDac(const ValSig* big, 
                             const ValSig* small,
                             const std::vector<ValSig>* bigSmallRatioN,
                             const std::vector<ValSig>* bigSmallRatioP) :
    m_bigSmallRatioN(0), m_bigSmallRatioP(0) 
  {
    if (big) {
      m_big = *big;
    }    else {
      m_big.setUndefined();
    }
    if (small) {
      m_small = *small;
    }    else {
      m_small.setUndefined();
    }
    if (bigSmallRatioN) {
      m_bigSmallRatioN = new std::vector<ValSig>(*bigSmallRatioN);
    }
    if (bigSmallRatioP) {
      m_bigSmallRatioP = new std::vector<ValSig>(*bigSmallRatioP);
    }
  }

  CalMevPerDac::~CalMevPerDac() {
    if (m_bigSmallRatioN) {
      delete m_bigSmallRatioN;
    }
    if (m_bigSmallRatioP) {
      delete m_bigSmallRatioP;
    }
  }

  const std::vector<ValSig>* CalMevPerDac::getBigSmallRatio(int face) const 
  {
    switch (face) {
    case idents::CalXtalId::POS:
      return m_bigSmallRatioP;
    case idents::CalXtalId::NEG:
      return m_bigSmallRatioN;
    default:
      return 0;
    }
  }


  void CalMevPerDac::update(RangeBase* other) {
    CalMevPerDac* otherConsts = dynamic_cast<CalMevPerDac* > (other);
    if (otherConsts) {
      m_big = otherConsts->m_big;
      m_small = otherConsts->m_small;
      if (m_bigSmallRatioN) {
        delete m_bigSmallRatioN;
        m_bigSmallRatioN = 0;
      }
      if (m_bigSmallRatioP) {
        delete m_bigSmallRatioP;
        m_bigSmallRatioP = 0;
      }
      if (otherConsts->m_bigSmallRatioN) {
        m_bigSmallRatioN = 
          new std::vector<ValSig>(*otherConsts->m_bigSmallRatioN);
      }
      if (otherConsts->m_bigSmallRatioP) {
        m_bigSmallRatioP = 
          new std::vector<ValSig>(*otherConsts->m_bigSmallRatioP);
      }
    }
  }

  // CalMevPerDacCol implementation
  CalMevPerDacCol::CalMevPerDacCol(unsigned nTowerRow, unsigned nTowerCol,
                                   unsigned nLayer, unsigned nXtal,
                                   unsigned nXpos) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, 1, 1, 0, nXpos)  {
    unsigned size = m_finder->getSize();
    
    CalMevPerDac* pConsts = new CalMevPerDac[size];
    for (unsigned ix = 0; ix < size; ix++) {
      m_ranges[ix] = pConsts; 
      ++pConsts;
    }

    if (nXpos) {  // Only possible values should be 0 or 1
      m_xpos = new Xpos;
    }
  }

  CalMevPerDacCol::~CalMevPerDacCol() {
    CalMevPerDac* pConsts = dynamic_cast<CalMevPerDac*>(m_ranges[0]);
    if (pConsts) delete [] pConsts;
    if (m_xpos) delete m_xpos;
  }

  bool CalMevPerDacCol::putRange(idents::CalXtalId id, unsigned range, 
                                 unsigned face, RangeBase* data) {
    if (!dynamic_cast<CalMevPerDac*>(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }

  bool CalMevPerDacCol::putRange(unsigned towerRow, unsigned towerCol, 
                               unsigned layer, unsigned xtal, unsigned range,
                               unsigned face, RangeBase* data) {
    if (!dynamic_cast<CalMevPerDac*>(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                  range, face, data);
  }

  const CLID& CalMevPerDacCol::classID() {return CLID_Calib_CAL_MevPerDac;}
}
