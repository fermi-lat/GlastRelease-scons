// $Header$

/**
      @file AncFinder.cxx

      AncFinder implementation.  Class "knows" how to find anc. data
      for a particular channel,layer, module.

       Calibration data is kept in what is logically a 3-dimensional
       array,
                   data[iMod, iLay, iChan]
       (for tagger; qdc has no layer field, handled by 
        always setting iLay=0)
       
       but will be stored as a singly-dimensioned array or
       perhaps vector.  
       Data for NA (not-attached) is kept in a separate singly-dimensioned
       array
*/
#include "AncFinder.h"

namespace CalibData {
  AncFinder::AncFinder(unsigned nMod, unsigned nLay, unsigned nChan) :
    : m_mod(nMod), m_lay(nLay), m_chan(nChan)
  {
    m_c0 = m_chan;
    m_c1 = m_c0 * m_lay;
  }


  bool AncFinder::equals(const AncFinder& other) const {
    return
      ((m_mod    == other.m_mod)    &&
       (m_lay     == other.m_lay)     &&
       (m_chan     == other.m_chan)   );
  }
}

