// $Header$
/** 
      @file TkrFinder.cxx

      TkrFinder implementation.  Class "knows" how to find Tkr. data
      for a particular tower, tray, top/bot layer (and perhaps FE chip)

       Calibration data is kept in what is logically a 3-dimensional
       array,
                   data[iTower,iUniplane, iFechip]
       but will be stored as a singly-dimensioned array or
       perhaps vector.  
*/
#include "CalibData/Tkr/TkrFinder.h"

namespace CalibData {
  unsigned TkrFinder::findIx(const idents::TkrId& id, unsigned feChip) const {
    return findIx(id.getTowerY(), id.getTowerX(), id.getTray(), 
                  id.getBotTop(), feChip);
  }

  unsigned TkrFinder::findIx(unsigned towerRow, unsigned towerCol, 
                             unsigned tray, bool top, unsigned feChip)  const
  {
    unsigned iTower = m_towerCol*towerRow + towerCol;
    unsigned uni = 2*tray + (top ? 1 : 0);
    return feChip + m_c0*uni + m_c1*iTower;
  }

  bool TkrFinder::equals(const TkrFinder& other) const {
    return
      ((m_towerRow == other.m_towerRow) &&
       (m_towerCol == other.m_towerCol) &&
       (m_unilayer == other.m_unilayer) &&
       (m_feChip == other.m_feChip) );
  }

  void TkrFinder::init() {
    m_tower = m_towerRow * m_towerCol;
    if (m_feChip) {
      m_c0 = m_feChip;
    }
    else m_c0 = 1;
    m_c1 = m_c0 * m_unilayer;
  }
}
