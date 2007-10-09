// $Header$

#include "AcdFinder.h"
#include "CalibData/RangeBase.h" 
#include "CalibData/Acd/AcdCalibBase.h"
#include "CalibData/Acd/AcdCalibObj.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {

  const CLID AcdCalibBase::noCLID = 0;
  const CLID& AcdCalibBase::classID() {return noCLID;}
  AcdCalibBase::AcdCalibBase(const AcdCalibDescription& desc,
			     unsigned nFace, unsigned nRow, unsigned nCol, 
                             unsigned nNA, unsigned nPmt) : 
    m_finder(0),
    m_desc(&desc){
    cGuts(nFace, nRow, nCol, nNA, nPmt);
  }

  AcdCalibBase::~AcdCalibBase() {
    delete m_finder;
  }

  AcdCalibObj*  AcdCalibBase::get(idents::AcdId id, unsigned pmt) { 
    int ix = m_finder->findIx(id, pmt);
    if (ix >= 0) {  // real detector behind electronics
      unsigned uix = ix;
      if (uix < m_finder->getSize() ) { 
        return m_pmts[uix];
      }
    }
    else {
      unsigned uix = NAindex(ix);
      if (uix < m_finder->getNNASize()) {
        return m_NAs[uix];
      }
    }
    return 0;    // failed
  }


  bool AcdCalibBase::put(idents::AcdId id, unsigned pmt, 
			 AcdCalibObj& pmtCalib) {
    int ix = m_finder->findIx(id, pmt);
    AcdCalibObj* pDest;
    if (ix >= 0) {
      unsigned uix = ix;
      if (uix >= m_finder->getSize() ) return false;
      pDest = m_pmts[uix];
      delete pDest;
      m_pmts[uix] = &pmtCalib;
    }  else {
      unsigned uix = NAindex(ix);
      if (uix >= m_finder->getNNASize()) return false;
      pDest = m_NAs[uix];
      delete pDest;
      m_NAs[uix] = &pmtCalib;
    }
    return true;
  }

  StatusCode AcdCalibBase::update(CalibBase& other, MsgStream* log) {
    AcdCalibBase& other1 = dynamic_cast<AcdCalibBase& >(other);

    // Check that new data is dimensioned the same as old
    if (!(m_finder->equals(*(other1.m_finder)) ) ) {  // tilt!  
      (*log) << MSG::ERROR 
             << "AcdCalibBase::update failure: dimensioning unequal" << endreq;
      return StatusCode::FAILURE;
    }

    // Update CalibBase stuff
    StatusCode sc = CalibBase::update(other, log);
    if (sc != StatusCode::SUCCESS) return sc;

    unsigned n = m_finder->getSize();
    for (unsigned i = 0; i < n; i++) {
      AcdCalibObj* dest = m_pmts[i];
      dest->update(other1.m_pmts[i]);
    }

    n = m_finder->getNNASize();
    for (unsigned j = 0; j < n; j++) {
      AcdCalibObj* dest = m_NAs[j];
      dest->update(other1.m_NAs[j]);
    }
    return StatusCode::SUCCESS;
  }

  void AcdCalibBase::cGuts(unsigned nFace, unsigned nRow, unsigned nCol, 
                           unsigned nNA, unsigned nPmt) {

    m_finder = new AcdFinder(nFace, nRow, nCol, nNA, nPmt);
    unsigned n = m_finder->getSize();

    //    m_pR = new vector<RangeBase*>(n, 0);
    m_pmts.resize(n,0);
    m_NAs.resize(nNA * nPmt,0);
  }




}
