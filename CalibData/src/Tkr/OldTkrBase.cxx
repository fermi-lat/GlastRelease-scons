// $Header$

#include "CalibData/Tkr/OldTkrFinder.h"
#include "CalibData/RangeBase.h" 
#include "CalibData/Tkr/OldTkrBase.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {

  const CLID OldTkrBase::noCLID = 0;
  const CLID& OldTkrBase::classID() {return noCLID;}
  OldTkrBase::OldTkrBase(unsigned nTowerRow, unsigned nTowerCol, 
                   unsigned nTray, bool indirect)
    : m_finder(0), m_indirect(indirect), m_ixValid(false) {
    cGuts(nTowerRow, nTowerCol, nTray);
  }

  OldTkrBase::~OldTkrBase() {
    delete m_finder;
  }

  RangeBase*  OldTkrBase::getChannel(const idents::TkrId& id) {
    m_ixValid = false;
    m_ix = m_finder->findIx(id);
    if (m_ix < m_finder->getSize() ) { 
      m_ixValid = true;
      if (m_indirect) return m_ranges[m_ix];
      else return 0;
    }
    return 0;
  }
  RangeBase*  OldTkrBase::getChannel(unsigned towerRow, unsigned towerCol, 
                                  unsigned tray, bool top) { 
    idents::TkrId  id(towerCol, towerRow, tray, top);
    return getChannel(id);
  }

  bool OldTkrBase::putChannel(RangeBase* data, unsigned towerRow, 
                           unsigned towerCol, unsigned tray, bool top) {
    return putChannel(data, idents::TkrId(towerCol, towerRow, tray, top));

  }
  bool OldTkrBase::putChannel(RangeBase* data, const idents::TkrId& id) {
    m_ixValid = false;
    m_ix = m_finder->findIx(id);
    if (m_ix >= m_finder->getSize()) return false;
    m_ixValid= true;

    //  if not indirect, derived class has to do the rest
    if (!m_indirect) return true;  

    
    RangeBase* pDest = m_ranges[m_ix];
    if (!pDest) {  // call back derived class to make a place for itself
      data->makeNew(&m_ranges[m_ix]);
      pDest = m_ranges[m_ix];
    }

    pDest->update(data);
    return true;
  }

  StatusCode OldTkrBase::update(CalibBase& other, MsgStream* log) {
    OldTkrBase& other1 = dynamic_cast<OldTkrBase& >(other);

    unsigned n = m_finder->getSize();

    // Check that new data is dimensioned the same as old
    if (!(m_finder->equals(*(other1.m_finder)) ) ) {  // tilt!  
      (*log) << MSG::ERROR 
             << "OldTkrBase::update failure: dimensioning unequal" << endreq;
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

  void OldTkrBase::cGuts(unsigned nTowerRow, unsigned nTowerCol, 
                      unsigned nTray) {

    m_finder = new OldTkrFinder(nTowerRow, nTowerCol, nTray);
    unsigned n = m_finder->getSize();

    // had this before
    //    if (m_indirect) m_ranges.reserve(n);

    // Following should allocate and also initialize all entries to 0
    if (m_indirect) m_ranges.resize(n, 0);
  }


}
