//$Header$
#ifndef CalibData_TkrFinder_h
#define CalibData_TkrFinder_h
#include "idents/TkrId.h"

/** @class TkrFinder

   Utility class which knows how to find the index of the
   correct set of TKR calibration constants given id.

   Calibration data is kept in what is logically a 2-dimensional
   array,
   data[iTower,iUniplane]

   Within a unilayer, data may be organized by strip or by fe chip.
   There is also a small amount of per-tower information
   
   The base class keeps an array of 16 pointers to TkrTower
   objects (since there can never be more than 16), initialized
   to 0.  TkrTowers are allocated as needed.  
   
   Note convention for tower numbering follows
   (Ritz, Nordby) LAT-SS-00035-03-D1.  For the standard 4x4 array,
   this yields

                     12   13   14   15  | r3
                                        |
                      8    9   10   11  | r2
                                        |           ^
                      4    5    6    7  | r1        |
                                        |           |
                      0    1    2    3  | r0        
                                        |           Y
                     -------------------+
                     c0   c1   c2   c3              X -->

                            

   @author J. Bogart

*/

namespace CalibData {
  class TkrFinder {
  public: 
    /**
       Constructor. Input defines 
     */
    TkrFinder(unsigned nTowerRow, unsigned nTowerCol, unsigned nTray = 19,
              unsigned nFeChip = 0) : m_towerRow(nTowerRow), 
                                      m_towerCol(nTowerCol), 
                                      m_nUni(2*nTray)
    {m_nTower = m_towerRow * m_towerCol;}

    ~TkrFinder() {}

    unsigned findTowerIx(unsigned towerRow, unsigned towerCol)     { 
      return towerCol + m_towerCol*towerRow;              }

    unsigned findTowerIx(const idents::TkrId& id) {
      return findTowerIx(id.getTowerY(), id.getTowerX());}

    // Check that fields of supplied id are within our range
    bool okId(const idents::TkrId& id) {
      return 
        (   (id.getTowerX() < m_towerCol) &&
            (id.getTowerY() < m_towerRow) &&
            (2*id.getTray() <= m_nUni)           );
    }

    unsigned findUniIx(const idents::TkrId& id)                            {
      return (2*id.getTray() + (id.getBotTop() ? 1 : 0) ); }

    unsigned findTowerRow(unsigned towerIx) {return towerIx / m_towerCol;}

    unsigned findTowerCol(unsigned towerIx) {return towerIx % m_towerCol;}
    unsigned getNTowerRow() const {return m_towerRow;}
    unsigned getNTowerCol() const {return m_towerCol;}
    unsigned getNUnilayer() const {return m_nUni;}

    bool equals(const TkrFinder& other) const {
    return
      ((m_towerRow == other.m_towerRow) &&
       (m_towerCol == other.m_towerCol) &&
       (m_nUni == other.m_nUni) );
  }

  private:
    TkrFinder() : m_towerRow(4), m_towerCol(4), m_nUni(38)
                   {m_nTower=m_towerRow * m_towerCol;}  // uses all defaults
    unsigned m_towerRow;
    unsigned m_towerCol;
    unsigned m_nTower;
    unsigned m_nUni;

  };
}    

#endif
