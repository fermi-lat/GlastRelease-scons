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
  CalFinder::CalFinder(unsigned nTower, unsigned nLayer, unsigned nColumn, 
                       unsigned nFace=2, unsigned nRange=4) : 
    m_tower(nTower), m_layer(nLayer),
    m_column(nColumn), m_face(nFace),
    m_range(nRange) {
    // Compute constants needed to find our way quickly in the array.
    m_c0 = m_face;
    m_c1 = m_c0 * m_range;
    m_c2 = m_c1 * m_column;
    m_c3 = m_c2 * m_layer;
  }

}

