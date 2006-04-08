// $Header$

/**
      @file AcdFinder.cxx

      AcdFinder implementation.  Class "knows" how to find acd. data
      for a particular tile-or-ribbon, pmt, range.

       Calibration data for real detectors (tiles and ribbons) is kept 
       in what is logically a 4-dimensional
       array,
                   data[iFace, iRow, iCol, iPmt]   for tile
                   data[iOrient, 0, iRibbon, iPmt] for ribbon
       
       but will be stored as a singly-dimensioned array or
       perhaps vector.  
       Data for NA (not-attached) is kept in a separate singly-dimensioned
       array
*/
#include "AcdFinder.h"

namespace CalibData {
  AcdFinder::AcdFinder(unsigned nFace, unsigned nRow, unsigned nCol, 
                       unsigned nNA, unsigned nPmt /*, unsigned  nDacCol */)
    : m_face(nFace), m_row(nRow), m_col(nCol), m_NA(nNA), m_pmt(nPmt)
    /* ,m_dacCol(nDacCol) */ {
    // Compute constants needed to find our way quickly in the array.
    m_c0 = m_pmt;
    m_c1 = m_c0 * m_col;
    m_c2 = m_c1 * m_row;
  }

  int AcdFinder::findIx(idents::AcdId id, unsigned pmt) const {

    if (id.na() == 1) {
      return pmt + m_pmt * NAencoded(id.colLike());
    }

    unsigned face = id.faceLike();
    unsigned row = id.rowLike();
    unsigned col = id.colLike();

    return 
      findIx(face, row, col, pmt);
    }

  bool AcdFinder::equals(const AcdFinder& other) const {
    return
      ((m_face    == other.m_face)    &&
       (m_row     == other.m_row)     &&
       (m_col     == other.m_col)     &&
       (m_NA      == other.m_NA)     &&
       (m_pmt     == other.m_pmt)   );
  }
}

