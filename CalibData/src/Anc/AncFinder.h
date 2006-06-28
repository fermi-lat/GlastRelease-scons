//$Header$
#ifndef CalibData_AncFinder_h
#define CalibData_AncFinder_h

/** @class AncFinder

   Utility class which knows how to find the index of the
   correct set of Anc calibration constants given module, layer, chan

   @author J. Bogart

*/

namespace CalibData {

  class AncFinder {
  public: 
    AncFinder(unsigned nMod, unsigned nLay, unsigned nChan);

    ~AncFinder() {}

    /**
       Find index of "regular" 
     */
    int findIx(unsigned face, unsigned row, unsigned col, 
                    unsigned pmt) const {
      return chan + m_c0*lay + m_c1*mod;
    }

    unsigned getSize() const {return m_c1*m_mod;}

    bool equals(const AncFinder& other) const;

  private:
    unsigned m_mod;
    unsigned m_lay;
    unsigned m_chan;

    unsigned m_c0;
    unsigned m_c1;
  };
}    



#endif
