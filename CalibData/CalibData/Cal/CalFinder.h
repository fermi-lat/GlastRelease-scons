//$Header$
#ifndef CalibData_CalFinder_h
#define CalibData_CalFinder_h

/** @class CalFinder

   Utility class which knows how to find the index of the
   correct set of CAL calibration constants given id, range
   and (optional) face

   @author J. Bogart

*/
#include "idents/CalXtalId.h"

namespace CalibData {
  class CalFinder {
  public: 
    CalFinder(unsigned nTower, unsigned nLayer, unsigned nColumn, 
              unsigned nFace=2, unsigned nRange=4);


    ~CalFinder() {}

    unsigned findIx(unsigned tower, unsigned layer, unsigned column, 
                    unsigned range, unsigned face=0) const {
      return face + m_c0*range + m_c1*column + m_c2*layer + m_c3*tower;
    }

    unsigned findIx(idents::CalXtalId id, unsigned range, unsigned face=0) 
      const
    {
      return 
        findIx(id.getTower(), id.getLayer(), id.getColumn(), range, face);
    }

    unsigned getSize() const {return m_c3*m_tower;}

  private:
    unsigned m_tower;
    unsigned m_layer;
    unsigned m_column;
    unsigned m_face;
    unsigned m_range;

    unsigned m_c0;
    unsigned m_c1;
    unsigned m_c2;
    unsigned m_c3;

  };
}    



#endif
