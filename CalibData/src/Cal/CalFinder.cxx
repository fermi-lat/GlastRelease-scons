// $Header$

/**
      @file CalFinder.cxx

      CalFinder implementation.  Class "knows" how to find cal. data
      for a particular crystal, range, face.

       Calibration data is kept in what is logically a 5-dimensional
       array,
                   data[iTower,iLayer,iCol, iRange, iFace]
       but will be stored as a singly-dimensioned array or
       perhaps vector.  
*/
#include "CalibData/Cal/CalFinder.h"

namespace CalibData {
  CalFinder::CalFinder(unsigned nTowerRow, unsigned nTowerCol, unsigned nLayer,
                       unsigned nXtal, unsigned nFace, unsigned nRange,
                       unsigned nDacCol, unsigned nXpos) : 
    m_towerRow(nTowerRow), m_towerCol(nTowerCol), 
    m_tower(nTowerRow*nTowerCol), m_layer(nLayer),
    m_xtal(nXtal), m_face(nFace), m_range(nRange),
    m_dacCol(nDacCol), m_xpos(nXpos) {
    // Compute constants needed to find our way quickly in the array.
    m_c0 = m_face;
    m_c1 = m_c0 * m_range;
    m_c2 = m_c1 * m_xtal;
    m_c3 = m_c2 * m_layer;
  }

  unsigned CalFinder::findIx(idents::CalXtalId id, unsigned range, 
                             unsigned face)       const
  {
    unsigned iTow = id.getTower();
    unsigned iRow = iTow / m_towerCol;
    unsigned iCol = iTow - iRow*m_towerCol;
    return 
      findIx(iRow, iCol, id.getLayer(), id.getColumn(), range, face);
    }

  bool CalFinder::checkIx(idents::CalXtalId id, unsigned range, 
                          unsigned face) {
    unsigned iTow = id.getTower();
    unsigned iRow = iTow / m_towerCol;
    unsigned iCol = iTow - iRow*m_towerCol;
    return 
      checkIx(iRow, iCol, id.getLayer(), id.getColumn(), range, face);
  }

  bool CalFinder::equals(const CalFinder& other) const {
    return
      ((m_towerRow == other.m_towerRow) &&
       (m_towerCol == other.m_towerCol) &&
       (m_layer    == other.m_layer)    &&
       (m_xtal     == other.m_xtal)     &&
       (m_face     == other.m_face)     &&
       (m_range    == other.m_range)    &&
       (m_dacCol   == other.m_dacCol)   &&
       (m_xpos     == other.m_xpos)  );
  }
}

