// $Header$

#include "CalibData/Acd/AcdFinder.h"
#include "CalibData/RangeBase.h" 
#include "CalibData/Acd/AcdCalibBase.h"
// #include "CalibData/DacCol.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {

  const CLID AcdCalibBase::noCLID = 0;
  const CLID& AcdCalibBase::classID() {return noCLID;}
  AcdCalibBase::AcdCalibBase(unsigned nFace, unsigned nRow, unsigned nCol, 
                             unsigned nPmt, unsigned nRange, 
                             unsigned /* nDacCol */ ) : m_finder(0) {

    cGuts(nFace, nRow, nCol, nPmt, nRange);
  }

  AcdCalibBase::~AcdCalibBase() {
    delete m_finder;
  }

  RangeBase*  AcdCalibBase::getRange(idents::AcdId id, unsigned pmt, 
                                     unsigned range) {
    unsigned ix = m_finder->findIx(id, pmt, range);
    if (ix < m_finder->getSize() ) { 
      return m_ranges[ix];
    }
    else return 0;
  }

  bool AcdCalibBase::putRange(unsigned face, unsigned row, unsigned col, 
                              unsigned pmt, unsigned range,
                              RangeBase* data) {
    unsigned ix = m_finder->findIx(face, row, col, pmt, range);
    if (ix >= m_finder->getSize() ) return false;

    RangeBase* pDest = m_ranges[ix];

    pDest->update(data);
    return true;
  }

  bool AcdCalibBase::putRange(idents::AcdId id, unsigned pmt, unsigned range, 
                              RangeBase* data) {
    unsigned ix = m_finder->findIx(id, pmt, range);
    if (ix >= m_finder->getSize() ) return false;

    RangeBase* pDest = m_ranges[ix];

    pDest->update(data);
    return true;
  }
  /*
  bool CalCalibBase::putDacCol(unsigned range, DacCol* dacs) {
    if (range >= m_dacCols.size()) return false;
    m_dacCols[range]->update(dacs);
    return true;
  }


  DacCol* CalCalibBase::getDacCol(unsigned range) {
    if (range >= m_dacCols.size()) return 0;
    return m_dacCols[range];
  }
  */

  StatusCode AcdCalibBase::update(CalibBase& other, MsgStream* log) {
    AcdCalibBase& other1 = dynamic_cast<AcdCalibBase& >(other);

    unsigned n = m_finder->getSize();

    // Check that new data is dimensioned the same as old
    if (!(m_finder->equals(*(other1.m_finder)) ) ) {  // tilt!  
      (*log) << MSG::ERROR 
             << "AcdCalibBase::update failure: dimensioning unequal" << endreq;
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

  void AcdCalibBase::cGuts(unsigned nFace, unsigned nRow, unsigned nCol, 
                           unsigned nPmt, unsigned nRange) {

    m_finder = new AcdFinder(nFace, nRow, nCol, nPmt, nRange);
    unsigned n = m_finder->getSize();

    //    m_pR = new vector<RangeBase*>(n, 0);
    m_ranges.reserve(n);

  }


}
