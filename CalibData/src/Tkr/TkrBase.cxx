// $Header$

#include "CalibData/Tkr/TkrFinder.h"
#include "CalibData/RangeBase.h" 
#include "CalibData/Tkr/TkrBase.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {

  const CLID TkrBase::noCLID = 0;
  const CLID& TkrBase::classID() {return noCLID;}
  TkrBase::TkrBase(unsigned nTowerRow, unsigned nTowerCol, 
                   unsigned nTray, unsigned nFeChip, bool indirect)
    : m_finder(0), m_indirect(indirect), m_ixValid(false) {
    cGuts(nTowerRow, nTowerCol, nTray, nFeChip);
  }

  TkrBase::~TkrBase() {
    delete m_finder;
  }

  RangeBase*  TkrBase::getChannel(const idents::TkrId& id, unsigned feChip) {
    m_ixValid = false;
    m_ix = m_finder->findIx(id, feChip);
    if (m_ix < m_finder->getSize() ) { 
      m_ixValid = true;
      if (m_indirect) return m_ranges[m_ix];
      else return 0;
    }
    return 0;
  }
  RangeBase*  TkrBase::getChannel(unsigned towerRow, unsigned towerCol, 
                                  unsigned tray, bool top, unsigned feChip) { 
    idents::TkrId  id(towerCol, towerRow, tray, top);
    return getChannel(id, feChip);
  }

  bool TkrBase::putChannel(RangeBase* data, unsigned towerRow, 
                           unsigned towerCol, unsigned tray, bool top,
                           unsigned feChip) {
    return putChannel(data, idents::TkrId(towerCol, towerRow, tray, top),
                      feChip);
  }
  bool TkrBase::putChannel(RangeBase* data, const idents::TkrId& id,
                           unsigned feChip) {
    m_ixValid = false;
    m_ix = m_finder->findIx(id, feChip);
    if (m_ix >= m_finder->getSize()) return false;
    m_ixValid= true;

    //  if not indirect, derived class has to do the rest
    if (!m_indirect) return true;  

    RangeBase* pDest = m_ranges[m_ix];
    pDest->update(data);
    return true;
  }

  StatusCode TkrBase::update(CalibBase& other, MsgStream* log) {
    TkrBase& other1 = dynamic_cast<TkrBase& >(other);

    unsigned n = m_finder->getSize();

    // Check that new data is dimensioned the same as old
    if (!(m_finder->equals(*(other1.m_finder)) ) ) {  // tilt!  
      (*log) << MSG::ERROR 
             << "TkrBase::update failure: dimensioning unequal" << endreq;
      return StatusCode::FAILURE;
    }

    // Update CalibBase stuff
    StatusCode sc = CalibBase::update(other, log);
    if (sc != StatusCode::SUCCESS) return sc;

    // Derived class has to do the work in case data is kept as vector of
    // values
    if (!m_indirect) return StatusCode::SUCCESS;

    unsigned i;
    for (i = 0; i < n; i++) {
      RangeBase* dest = m_ranges[i];
      dest->update(other1.m_ranges[i]);
    }
    return StatusCode::SUCCESS;
  }

  void TkrBase::cGuts(unsigned nTowerRow, unsigned nTowerCol, 
                      unsigned nTray, unsigned nChip) {

    m_finder = new TkrFinder(nTowerRow, nTowerCol, nTray, nChip);
    unsigned n = m_finder->getSize();

    //    m_pR = new vector<RangeBase*>(n, 0);
    if (m_indirect) m_ranges.reserve(n);
  }


}
