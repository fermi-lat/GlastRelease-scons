// $Header$

/**
  @file  CalTholdMuon.cxx

  Implementation for calorimeter muon threshold calibration, 
  organized by crystal face.  Some of the components are per xtal face, others
  are per range, per face.

  Implementation for two classes may be found here.  @a CalTholdMuon is a class
  for per-crystal-face information.  @a CalTholdMuonCol contains CalTholdMuon 
  information for the entire detector.
*/
#include "CalibData/Cal/CalTholdMuon.h"

namespace CalibData {
  // CalTholdMuon implementation
  CalTholdMuon::CalTholdMuon(const ValSig* FLE, 
                             const ValSig* FHE,
                             const std::vector<ValSig>* ped) : m_ped(0)
  {
    if (FLE) {
      m_FLE = *FLE;
    }    else {
      m_FLE.setUndefined();
    }
    if (FHE) {
      m_FHE = *FHE;
    }    else {
      m_FHE.setUndefined();
    }
    if (ped) {
      m_ped = new std::vector<ValSig>(*ped);
    }
  }

  CalTholdMuon::~CalTholdMuon() {
    if (m_ped) {
      delete m_ped;
    }
  }

  const ValSig* CalTholdMuon::getPed(int range) const 
  {
    switch (range) {
    case idents::CalXtalId::LEX8:
    case idents::CalXtalId::LEX1:
    case idents::CalXtalId::HEX8:
    case idents::CalXtalId::HEX1:
      return &(*m_ped)[range];

    default:
      return 0;
    }
  }


  void CalTholdMuon::update(RangeBase* other) {
    CalTholdMuon* otherConsts = dynamic_cast<CalTholdMuon* > (other);
    if (otherConsts) {
      m_FLE = otherConsts->m_FLE;
      m_FHE = otherConsts->m_FHE;

      if (m_ped) {
        delete m_ped;
        m_ped = 0;
      }
      if (otherConsts->m_ped) {
        m_ped = new std::vector<ValSig>(*otherConsts->m_ped);
      }
    }
  }

  // CalTholdMuonCol implementation
  CalTholdMuonCol::CalTholdMuonCol(unsigned nTowerRow, unsigned nTowerCol,
                                   unsigned nLayer, unsigned nXtal, 
                                   unsigned nFace) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, nFace, 1)  {
    unsigned size = m_finder->getSize();
    
    CalTholdMuon* pConsts = new CalTholdMuon[size];
    for (unsigned ix = 0; ix < size; ix++) {
      m_ranges[ix] = pConsts; 
      ++pConsts;
    }
  }

  CalTholdMuonCol::~CalTholdMuonCol() {
    CalTholdMuon* pConsts = dynamic_cast<CalTholdMuon*>(m_ranges[0]);
    if (pConsts) delete [] pConsts;
  }

  bool CalTholdMuonCol::putRange(idents::CalXtalId id, unsigned range, 
                                 unsigned face, RangeBase* data) {
    if (!dynamic_cast<CalTholdMuon*>(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }

  bool CalTholdMuonCol::putRange(unsigned towerRow, unsigned towerCol, 
                                 unsigned layer, unsigned xtal, unsigned range,
                                 unsigned face, RangeBase* data) {
    if (!dynamic_cast<CalTholdMuon*>(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                  range, face, data);
  }

  const CLID& CalTholdMuonCol::classID() {return CLID_Calib_CAL_TholdMuon;}
}
