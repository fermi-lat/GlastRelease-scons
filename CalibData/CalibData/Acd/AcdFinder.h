//$Header$
#ifndef CalibData_AcdFinder_h
#define CalibData_AcdFinder_h

/** @class AcdFinder

   Utility class which knows how to find the index of the
   correct set of Acd calibration constants given id, pmt and range.
   Acd id's follow standard conventions.
                            

   @author J. Bogart

*/
#include "idents/AcdId.h"

namespace CalibData {
  class AcdFinder {
  public: 
    AcdFinder(unsigned nFace, unsigned nRow, unsigned nCol,
              unsigned nPmt=2, unsigned nRange=2,
              unsigned nDacCol=0);


    ~AcdFinder() {}

    unsigned findIx(unsigned face, unsigned row, unsigned col, 
                    unsigned pmt, unsigned range) const {
      //      unsigned iTower = m_towerCol*towerRow + towerCol;
      //      return face + m_c0*range + m_c1*xtal + m_c2*layer + m_c3*iTower;
      return range + m_c0*pmt + m_c1*col + m_c2*row + m_c3*face;
    }

    unsigned findIx(idents::AcdId id, unsigned pmt, unsigned range)
      const;

    unsigned getSize() const {return m_c3*m_face;}

    unsigned getNRange() const {return m_range;}

    //    unsigned getNDacCol() const {return m_dacCol;}

    bool equals(const AcdFinder& other) const;

  private:
    unsigned m_face;
    unsigned m_row;
    unsigned m_col;
    unsigned m_pmt;
    unsigned m_range;
    // unsigned m_dacCol;

    unsigned m_c0;
    unsigned m_c1;
    unsigned m_c2;
    unsigned m_c3;

  };
}    



#endif
