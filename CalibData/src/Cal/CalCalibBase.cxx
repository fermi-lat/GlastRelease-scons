// $Header$

#include "CalibData/Cal/CalFinder.h"
#include "CalibData/Cal/RangeBase.h" 
#include "CalibData/Cal/CalCalibBase.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {

  const CLID CalCalibBase::noCLID = 0;
  const CLID& CalCalibBase::classID() {return noCLID;}
  CalCalibBase::CalCalibBase(unsigned nTowerRow, unsigned nTowerCol, 
                             unsigned nLayer, unsigned nXtal, unsigned nFace, 
                             unsigned nRange) : m_finder(0) {
    m_finder = new CalFinder(nTowerRow, nTowerCol, nLayer, nXtal, nFace, 
                             nRange);
    unsigned n = m_finder->getSize();

    //    m_pR = new vector<RangeBase*>(n, 0);
    m_ranges.reserve(n);
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

  StatusCode CalCalibBase::update(CalibBase& other, MsgStream* log) {
    CalCalibBase& other1 = dynamic_cast<CalCalibBase& >(other);

    unsigned n = m_finder->getSize();

    // Make a simple but insufficient check that the new data is
    // structured like the old
    if (n != other1.m_finder->getSize() ) {  // tilt!  
      (*log) << MSG::ERROR 
             << "CalCalibBase::update failure: sizes unequal" << endreq;
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

}
