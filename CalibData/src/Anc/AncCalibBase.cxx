// $Header$

#include "AncFinder.h"
#include "CalibData/RangeBase.h" 
#include "CalibData/Anc/AncCalibBase.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {

  const CLID AncCalibBase::noCLID = 0;
  const CLID& AncCalibBase::classID() {return noCLID;}
  AncCalibBase::AncCalibBase(unsigned nMod, nLay, nChan) : m_finder(0) {
    cGuts(nMod, nLay, nChan);
  }

  AncCalibBase::~AncCalibBase() {
    delete m_finder;
  }

  RangeBase*  AncCalibBase::getChan(unsigned iMod, unsigned iLayer,
                                    unsigned iChan) {
    int ix = m_finder->findIx(iMod, iLayer, iChan);
    unsigned uix = ix;
    if (uix < m_finder->getSize() ) { 
      return m_chans[uix];
    }
    return 0;    // failed
  }

  /* Do we need this?  */
  bool AncCalibBase::putChan(unsigned mod, unsigned layer, unsigned chan,
                             RangeBase* data) {
    int ix = m_finder->findIx(mod, layer, chan);
    unsigned uix = ix;
    if (uix >= m_finder->getSize() ) return false;

    RangeBase* pDest = m_chans[uix];

    pDest->update(data);
    return true;
  }


  StatusCode AncCalibBase::update(CalibBase& other, MsgStream* log) {
    AncCalibBase& other1 = dynamic_cast<AncCalibBase& >(other);


    // Check that new data is dimensioned the same as old
    if (!(m_finder->equals(*(other1.m_finder)) ) ) {  // tilt!  
      (*log) << MSG::ERROR 
             << "AncCalibBase::update failure: dimensioning unequal" << endreq;
      return StatusCode::FAILURE;
    }

    // Update CalibBase stuff
    StatusCode sc = CalibBase::update(other, log);
    if (sc != StatusCode::SUCCESS) return sc;

    unsigned n = m_finder->getSize();
    for (unsigned i = 0; i < n; i++) {
      RangeBase* dest = m_chans[i];
      dest->update(other1.m_chans[i]);
    }
    return StatusCode::SUCCESS;
  }

  void AncCalibBase::cGuts(unsigned nMod, unsigned nLay, unsigned nChan) {
    m_finder = new AncFinder(nMod, nLay, nChan);
    unsigned n = m_finder->getSize();
    m_chans.reserve(n);
  }

}
