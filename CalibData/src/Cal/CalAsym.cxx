// $Header$

/**
  @file  CalAsym.cxx

  Implementation for new calorimeter asymmetry calibration, consisting of 
  asymmetry constants for all 4 combinations of diodes, one from each face.
  For each such combination store an array of constants corresponding to
  a (fixed for all crystals) array of x positions along the crystal.
  The array of positions is kept as part of the calibration.

  Implementation for two classes may be found here.  @a CalAsym is a class
  for per-crystal information.  @a CalAsymCol contains CalAsym 
  information for the entire detector and also the array of x positions.
*/
#include "CalibData/Cal/CalAsym.h"
#include "CalibData/Cal/Xpos.h"

namespace CalibData {
  // CalAsym implementation
  CalAsym::CalAsym(const std::vector<ValSig>* big,
                   const std::vector<ValSig>* small,
                   const std::vector<ValSig>* nSmallPBig,
                   const std::vector<ValSig>* pSmallNBig) :
    m_big(0), m_small(0), m_nSmallPBig(0), m_pSmallNBig(0)
  {
    if (big) {
      m_big = new std::vector<ValSig>(*big);
    }   

    if (small) {
      m_small = new std::vector<ValSig>(*small);
    } 

    if (nSmallPBig) {
      m_nSmallPBig = new std::vector<ValSig>(*nSmallPBig);
    }
    if (pSmallNBig) {
      m_pSmallNBig = new std::vector<ValSig>(*pSmallNBig);
    }
  }

  CalAsym::~CalAsym() {
    if (m_big) delete m_big;
    if (m_small) delete m_small;
    if (m_nSmallPBig) delete m_nSmallPBig;
    if (m_pSmallNBig) delete m_pSmallNBig;
  }

  void CalAsym::update(RangeBase* other) {
    CalAsym* otherConsts = dynamic_cast<CalAsym* > (other);
    if (otherConsts) {
      if (m_big) {
        delete m_big;
        m_big = 0;
      }
      if (m_small) {
        delete m_small;
        m_small = 0;
      }
      if (m_nSmallPBig) {
        delete m_nSmallPBig;
        m_nSmallPBig = 0;
      }
      if (m_pSmallNBig) {
        delete m_pSmallNBig;
        m_pSmallNBig = 0;
      }

      if (otherConsts->m_big) {
        m_big = new std::vector<ValSig>(*otherConsts->m_big);
      }
      if (otherConsts->m_small) {
        m_small = new std::vector<ValSig>(*otherConsts->m_small);
      }

      if (otherConsts->m_nSmallPBig) {
        m_nSmallPBig = 
          new std::vector<ValSig>(*otherConsts->m_nSmallPBig );
      }
      if (otherConsts->m_pSmallNBig) {
        m_pSmallNBig = 
          new std::vector<ValSig>(*otherConsts->m_pSmallNBig);
      }
    }
  }

  // CalAsymCol implementation
  CalAsymCol::CalAsymCol(unsigned nTowerRow, unsigned nTowerCol,
                                   unsigned nLayer, unsigned nXtal,
                                   unsigned nXpos) :
    CalCalibBase(nTowerRow, nTowerCol, nLayer, nXtal, 1, 1, 0, nXpos)  {
    unsigned size = m_finder->getSize();
    
    CalAsym* pConsts = new CalAsym[size];
    for (unsigned ix = 0; ix < size; ix++) {
      m_ranges[ix] = pConsts; 
      ++pConsts;
    }

    if (nXpos) {  // Only possible values should be 0 or 1
      m_xpos = new Xpos;
    }
  }

  CalAsymCol::~CalAsymCol() {
    CalAsym* pConsts = dynamic_cast<CalAsym*>(m_ranges[0]);
    if (pConsts) delete [] pConsts;
    if (m_xpos) delete m_xpos;
  }

  bool CalAsymCol::putRange(idents::CalXtalId id, unsigned range, 
                            unsigned face, RangeBase* data) {
    if (!dynamic_cast<CalAsym*>(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(id, range, face, data);
  }

  bool CalAsymCol::putRange(unsigned towerRow, unsigned towerCol, 
                            unsigned layer, unsigned xtal, unsigned range,
                            unsigned face, RangeBase* data) {
    if (!dynamic_cast<CalAsym*>(data)) return false;

    // Otherwise go ahead and let base class handle it
    return CalCalibBase::putRange(towerRow, towerCol, layer, xtal,
                                  range, face, data);
  }

  const CLID& CalAsymCol::classID() {return CLID_Calib_CAL_Asym;}
}
