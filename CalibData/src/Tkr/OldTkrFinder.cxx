// $Header$
/** 
      @file OldTkrFinder.cxx

      OldTkrFinder implementation.  Class "knows" how to find Tkr. data
      for a particular tower, tray, top/bot layer (and perhaps FE chip)

       Calibration data is kept in what is logically a 3-dimensional
       array,
                   data[iTower,iUniplane, iFechip]
       but will be stored as a singly-dimensioned array or
       perhaps vector.  
*/
#include "CalibData/Tkr/OldTkrFinder.h"

namespace CalibData {
  unsigned OldTkrFinder::findIx(const idents::TkrId& id) const {
    return findIx(id.getTowerY(), id.getTowerX(), id.getTray(), 
                  id.getBotTop());
  }

  unsigned OldTkrFinder::findIx(unsigned towerRow, unsigned towerCol, 
                             unsigned tray, bool top)  const
  {
    unsigned iTower = m_towerCol*towerRow + towerCol;
    unsigned uni = 2*tray + (top ? 1 : 0);
    return uni + m_unilayer*iTower;
  }

  bool OldTkrFinder::equals(const OldTkrFinder& other) const {
    return
      ((m_towerRow == other.m_towerRow) &&
       (m_towerCol == other.m_towerCol) &&
       (m_unilayer == other.m_unilayer) );
  }

  void OldTkrFinder::init() {
    m_tower = m_towerRow * m_towerCol;
    //    m_c0 = 1;
    //    m_c1 = m_c0 * m_unilayer;
  }
}
