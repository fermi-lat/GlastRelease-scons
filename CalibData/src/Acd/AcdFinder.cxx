// $Header$

/**
      @file AcdFinder.cxx

      AcdFinder implementation.  Class "knows" how to find acd. data
      for a particular tile-or-ribbon, pmt, range.

       Calibration data is kept in what is logically a 5-dimensional
       array,
                   data[iFace, iRow, iCol, iPmt, iRange]
       but will be stored as a singly-dimensioned array or
       perhaps vector.  
*/
#include "CalibData/Acd/AcdFinder.h"

namespace CalibData {
  AcdFinder::AcdFinder(unsigned nFace, unsigned nRow, unsigned nCol, 
                       unsigned nPmt, unsigned nRange,
                       unsigned /* nDacCol */) : 
    m_face(nFace), m_row(nRow), m_col(nCol), 
    m_pmt(nPmt), m_range(nRange)
    /* ,m_dacCol(nDacCol) */ {
    // Compute constants needed to find our way quickly in the array.
    m_c0 = m_range;
    m_c1 = m_c0 * m_pmt;
    m_c2 = m_c1 * m_col;
    m_c3 = m_c2 * m_row;
  }

  unsigned AcdFinder::findIx(idents::AcdId id, unsigned pmt, 
                             unsigned range) const { 

    unsigned face = id.face();
    unsigned row=0;
    unsigned col;
    if (id.tile() ) {
      row = id.row();
      col = id.column();
    }
    else {
      col = id.ribbonNum();
    }

    return 
      findIx(face, row, col, pmt, range);
    }

  bool AcdFinder::equals(const AcdFinder& other) const {
    return
      ((m_face    == other.m_face)    &&
       (m_row     == other.m_row)     &&
       (m_col     == other.m_col)     &&
       (m_pmt     == other.m_pmt)     &&
       (m_range    == other.m_range) );
    /*  && (m_dacCol   == other.m_dacCol) ); */
  }
}

