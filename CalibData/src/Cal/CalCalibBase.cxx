// $Header$

#include "CalibData/Cal/CalFinder.h"
#include "CalibData/RangeBase.h" 
#include "CalibData/Cal/CalCalibBase.h"
#include "CalibData/DacCol.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {

  const CLID CalCalibBase::noCLID = 0;
  const CLID& CalCalibBase::classID() {return noCLID;}
  CalCalibBase::CalCalibBase(unsigned nTowerRow, unsigned nTowerCol, 
                             unsigned nLayer, unsigned nXtal, unsigned nFace, 
                             unsigned nRange, unsigned nDacCol) : m_finder(0) {

    cGuts(nTowerRow, nTowerCol, nLayer, nXtal, nFace, nRange, nDacCol);
  }

  CalCalibBase::~CalCalibBase() {
    delete m_finder;
  }

  RangeBase*  CalCalibBase::getRange(idents::CalXtalId id, unsigned range, 
                                     unsigned face) {
    unsigned ix = m_finder->findIx(id, range, face);
    if (ix < m_finder->getSize() ) { 
      return m_ranges[ix];
    }
    else return 0;
  }
  RangeBase*  CalCalibBase::getRange(unsigned towerRow, unsigned towerCol, 
                                     unsigned layer, unsigned xtal, 
                                     unsigned range, unsigned face) {
    unsigned ix = m_finder->findIx(towerRow, towerCol, layer, xtal, 
                                   range, face);
    if (ix < m_finder->getSize() ) { 
      return m_ranges[ix];
    }
    else return 0;
  }

  bool CalCalibBase::putRange(unsigned towerRow, unsigned towerCol, 
                              unsigned layer, unsigned xtal, unsigned range,
                              unsigned face, RangeBase* data) {
    unsigned ix = m_finder->findIx(towerRow, towerCol, layer, xtal, 
                                   range, face);
    if (ix >= m_finder->getSize() ) return false;

    RangeBase* pDest = m_ranges[ix];

    pDest->update(data);
    return true;
  }

  bool CalCalibBase::putRange(idents::CalXtalId id, unsigned range, 
                              unsigned face, RangeBase* data) {
    unsigned ix = m_finder->findIx(id, range, face);
    if (ix >= m_finder->getSize() ) return false;

    RangeBase* pDest = m_ranges[ix];

    pDest->update(data);
    return true;
  }

  bool CalCalibBase::putDacCol(unsigned range, DacCol* dacs) {
    if (range >= m_dacCols.size()) return false;
    m_dacCols[range]->update(dacs);
    return true;
  }


  DacCol* CalCalibBase::getDacCol(unsigned range) {
    if (range >= m_dacCols.size()) return 0;
    return m_dacCols[range];
  }


  StatusCode CalCalibBase::update(CalibBase& other, MsgStream* log) {
    CalCalibBase& other1 = dynamic_cast<CalCalibBase& >(other);

    unsigned n = m_finder->getSize();

    // Check that new data is dimensioned the same as old
    if (!(m_finder->equals(*(other1.m_finder)) ) ) {  // tilt!  
      (*log) << MSG::ERROR 
             << "CalCalibBase::update failure: dimensioning unequal" << endreq;
      return StatusCode::FAILURE;
    }

    // Update CalibBase stuff
    StatusCode sc = CalibBase::update(other, log);
    if (sc != StatusCode::SUCCESS) return sc;

    unsigned i;
    for (i = 0; i < n; i++) {
      RangeBase* dest = m_ranges[i];
      dest->update(other1.m_ranges[i]);
    }
    return StatusCode::SUCCESS;
  }

  void CalCalibBase::cGuts(unsigned nTowerRow, unsigned nTowerCol, 
                           unsigned nLayer, unsigned nXtal, 
                           unsigned nFace, unsigned nRange, 
                           unsigned nDacCol) {

    m_finder = new CalFinder(nTowerRow, nTowerCol, nLayer, nXtal, nFace, 
                             nRange, nDacCol);
    unsigned n = m_finder->getSize();

    //    m_pR = new vector<RangeBase*>(n, 0);
    m_ranges.reserve(n);

  }


}
