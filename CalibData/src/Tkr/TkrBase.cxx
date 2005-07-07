// $Header$

#include "CalibData/Tkr/UniBase.h" 
#include "CalibData/Tkr/TkrBase.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {

  const CLID TkrBase::noCLID = 0;
  const CLID& TkrBase::classID() {return noCLID;}
  TkrBase::TkrBase(unsigned nTowerRow, unsigned nTowerCol, 
                   unsigned nTray, bool indirect)
    : m_finder(0), m_factory(0), m_indirect(indirect)  {
    cGuts(nTowerRow, nTowerCol, nTray);
  }

  TkrBase::~TkrBase() {
    delete m_finder;
    for (unsigned ix = 0; ix < TKRBASE_MAXTOWER; ix++) {
      if (m_towers[ix]) delete m_towers[ix];
      m_towers[ix] = 0;
    }
  }

  TkrBase::TkrTower::~TkrTower() {
    unsigned nUni = m_unis.size();
    for (unsigned ix = 0; ix < nUni; ix++) {
      if (m_unis[ix]) {
        delete m_unis[ix];
        m_unis[ix] = 0;
      }
    }
  }

  UniBase*  TkrBase::getUni(const idents::TkrId& id) {
    if (!m_finder->okId(id)) return 0;

    unsigned iTower = m_finder->findTowerIx(id);
    unsigned iUni = m_finder->findUniIx(id);

    TkrTower* pTower = m_towers[iTower];
    if (!pTower) return 0;
    return pTower->m_unis[iUni];
  }

  UniBase*  TkrBase::getUni(unsigned towerRow, unsigned towerCol, 
                             unsigned tray, bool top)             { 
    idents::TkrId  id(towerCol, towerRow, tray, top);
    return getUni(id);
  }

  const std::string* TkrBase::getHwserial(unsigned towerRow, 
                                          unsigned towerCol) const {
    unsigned iTow=  m_finder->findTowerIx(towerRow, towerCol);

    if (m_towers[iTow]) return &(m_towers[iTow]->m_hwserial);
    else return 0;
  }


  bool TkrBase::putUni(UniBase* data, unsigned towerRow, 
                        unsigned towerCol, unsigned tray, bool top) {
    return putUni(data, idents::TkrId(towerCol, towerRow, tray, top));
  }
  bool TkrBase::putUni(UniBase* data, const idents::TkrId& id) {

    if (!m_finder->okId(id)) return false;
    unsigned towerIx = m_finder->findTowerIx(id);
    if (!m_towers[towerIx]) {
      // report error??  Should have already made tower
      // before trying to add uniplane
      return 0;
    }
    //  if not indirect, derived class has to do the rest
    if (!m_indirect) return true;  // might be able to dispense with m_indirect


    TkrTower* pTower = m_towers[towerIx];
    unsigned iUni = m_finder->findUniIx(id);
    UniBase* pDest = pTower->m_unis[iUni];

    if (!pDest) {  // call back derived class to make a place for itself
      pDest = m_factory->makeUni();

    }

    // Fill in data
    pDest->update(data);
    return true;
  }

  StatusCode TkrBase::update(CalibBase& other, MsgStream* log) {
    TkrBase& other1 = dynamic_cast<TkrBase& >(other);

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


    for (unsigned iTower = 0; iTower < TKRBASE_MAXTOWER; iTower++) {
      TkrTower* towerDest = m_towers[iTower];
      TkrTower* towerSrc  = other1.m_towers[iTower];
      if (towerSrc) {
        if (!towerDest) {
          m_towers[iTower] = towerDest 
            = makeTower(iTower, m_finder->getNUnilayer());
        }
        towerDest->m_iRow = towerSrc->m_iRow;
        towerDest->m_iCol = towerSrc->m_iCol;
        towerDest->m_hwserial = towerSrc->m_hwserial;
        unsigned nUni = towerSrc->m_unis.size();
        towerDest->resize(nUni);
        for (unsigned iUni=0; iUni < nUni; iUni++) {
          UniBase* pDest = towerDest->m_unis[iUni];
          if (!pDest) {
            towerDest->m_unis[iUni] = pDest = m_factory->makeUni();
          }
          pDest->update(towerSrc->m_unis[iUni]);
        }
      }
      else {   // no src, so shouldn't have dest
        if (towerDest) delete towerDest;
        m_towers[iTower] = 0;
      }
    }

    return StatusCode::SUCCESS;
  }

  unsigned TkrBase::getNTowerRow() const {return m_finder->getNTowerRow();}

  unsigned TkrBase::getNTowerCol() const {return m_finder->getNTowerCol();}

  unsigned TkrBase::getNUnilayer() const {return m_finder->getNUnilayer();}

  UniBase* UniFactoryBase::makeUni() {return new UniBase();}

  TkrBase::TkrTower* TkrBase::makeTower(unsigned iTow, unsigned nUni) {
    if (iTow < TKRBASE_MAXTOWER) {
      TkrTower* tow = m_towers[iTow] = new TkrTower();
      tow->resize(nUni);
      return tow;
    }
    return 0;
  }


  void TkrBase::cGuts(unsigned nTowerRow, unsigned nTowerCol, 
                      unsigned nTray) {

    m_finder = new TkrFinder(nTowerRow, nTowerCol, nTray);
    for (unsigned iTower = 0; iTower < TKRBASE_MAXTOWER; iTower++) {
      m_towers[iTower] = 0;
    }
    m_factory = new UniFactoryBase();

  }

  void TkrBase::TkrTower::resize(unsigned n)  {
    unsigned oldSize = m_unis.size();
    if (n != oldSize) {
      for (unsigned ix=0; ix < oldSize; ix++) {
        delete m_unis[ix];
      }
      m_unis.resize(n);
    }
  }


}
