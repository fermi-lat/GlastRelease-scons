// $Header$

/**
  @file  CalTholdCI.cxx

  Implementation for calorimeter charge injection threshold calibration, 
  organized by crystal face.  Some of the components are per face, others
  are per range, per face.

  Implementation for two classes may be found here.  @a CalTholdCI is a class
  for per-crystal-face information.  @a CalTholdCICol contains CalTholdCI 
  information for the entire detector.
*/
#include "CalibData/Cal/CalTholdCI.h"

namespace CalibData {
  // CalTholdCI implementation
  CalTholdCI::CalTholdCI(const std::vector<ValSig>* ULD,
                         const ValSig* FLE, 
                         const ValSig* FHE,
                         const ValSig* LAC,
                         const std::vector<ValSig>* ped) : m_ULD(0),
                                                           m_ped(0)
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
    if (LAC) {
      m_LAC = *LAC;
    }    else {
      m_LAC.setUndefined();
    }
    if (ULD) {
      m_ULD = new std::vector<ValSig>(*ULD);
    }
    if (ped) {
      m_ped = new std::vector<ValSig>(*ped);
    }
  }

  CalTholdCI::~CalTholdCI() {
    if (m_ped) delete m_ped;
    if (m_ULD) delete m_ULD;
  }

  const ValSig* CalTholdCI::getPed(int range) const 
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

  const ValSig* CalTholdCI::getULD(int range) const 
  {
    switch (range) {
    case idents::CalXtalId::LEX8:
    case idents::CalXtalId::LEX1:
    case idents::CalXtalId::HEX8:
    case idents::CalXtalId::HEX1:
      return &(*m_ULD)[range];

    default:
      return 0;
    }
  }

/*  Everything above has been adjusted for CalTholdCI */

/*       START HERE                    */

  void CalTholdCI::update(RangeBase* other) {
    CalTholdCI* otherConsts = dynamic_cast<CalTholdCI* > (other);
    if (otherConsts) {
      m_FLE = otherConsts->m_FLE;
      m_FHE = otherConsts->m_FHE;
      m_LAC = otherConsts->m_LAC;

      if (m_ULD) {
        delete m_ULD;
        m_ULD = 0;
      }
      if (m_ped) {
        delete m_ped;
        m_ped = 0;
      }
      if (otherConsts->m_ULD) {
        m_ULD = new std::vector<ValSig>(*otherConsts->m_ULD);
      }
      if (otherConsts->m_ped) {
        m_ped = new std::vector<ValSig>(*otherConsts->m_ped);
      }
    }
  }

  // CalTholdCICol implementation
  CalTholdCICol::CalTholdCICol(unsigned nTowerRow, unsigned nTowerCol,
                               unsigned nLayer, unsigned nXtal, 
                               unsigned nFace) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, nFace, 1)  {
    unsigned size = m_finder->getSize();
    
    CalTholdCI* pConsts = new CalTholdCI[size];
    for (unsigned ix = 0; ix < size; ix++) {
      m_ranges[ix] = pConsts; 
      ++pConsts;
    }
  }

  CalTholdCICol::~CalTholdCICol() {
    CalTholdCI* pConsts = dynamic_cast<CalTholdCI*>(m_ranges[0]);
    if (pConsts) delete [] pConsts;
  }

  bool CalTholdCICol::putRange(idents::CalXtalId id, unsigned range, 
                                 unsigned face, RangeBase* data) {
    if (!dynamic_cast<CalTholdCI*>(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }

  bool CalTholdCICol::putRange(unsigned towerRow, unsigned towerCol, 
                               unsigned layer, unsigned xtal, unsigned range,
                               unsigned face, RangeBase* data) {
    if (!dynamic_cast<CalTholdCI*>(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                  range, face, data);
  }

  const CLID& CalTholdCICol::classID() {return CLID_Calib_CAL_TholdCI;}
}
