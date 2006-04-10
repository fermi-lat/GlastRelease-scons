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

  // macros to translate to/from NA index.  (findIx returns an
  // encoded version so that same routine can handle NA index and
  // index for real PMTs
  // translate findIx return value to simple non-negative index
#define NAindex(ix) (-(ix) -1)
  // and back again.  Hey, it's the same transformation!
#define NAencoded(ix) (-(ix) -1)

  class AcdFinder {
  public: 
    AcdFinder(unsigned nFace, unsigned nRow, unsigned nCol, unsigned nNA,
              unsigned nPmt=2);
              //              unsigned nDacCol=0);

    ~AcdFinder() {}

    /**
       Find index of "regular" 
     */
    int findIx(unsigned face, unsigned row, unsigned col, 
                    unsigned pmt) const {
      return pmt + m_c0*col + m_c1*row + m_c2*face;
    }

    /**
        Returns index into "regular" PMT information or NA channels.
        If return is >= 0, refers to regular PMT.  If < 0, it's NA, 
        according to formula
            retValue = -((2*iNA + iPmt + 1))
        So index into NA array will be
            -retValue - 1
     */
    int findIx(idents::AcdId id, unsigned pmt) const;

    unsigned getSize() const {return m_c2*m_face;}
    unsigned getNNASize() const {return (m_NA*m_pmt);}

    bool equals(const AcdFinder& other) const;

  private:
    unsigned m_face;
    unsigned m_row;
    unsigned m_col;
    unsigned m_NA;
    unsigned m_pmt;
    //    unsigned m_range;
    // unsigned m_dacCol;

    unsigned m_c0;
    unsigned m_c1;
    unsigned m_c2;
    //    unsigned m_c3;

  };
}    



#endif
