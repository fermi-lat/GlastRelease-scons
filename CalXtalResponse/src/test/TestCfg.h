// $Header$
/** @file
    @author Z.Fewtrell

    
*/
#ifndef TestCfg_h
#define TestCfg_h

// LOCAL

// GLAST
#include "CalUtil/CalDefs.h"
#include "idents/CalXtalId.h"

// EXTLIB

// STD
#include <vector>
#include <set>

class ICalCalibSvc;

namespace CalXtalResponse {

  typedef std::set<CalUtil::TwrNum> TwrSet;

  /** class represents test parameters for single CalXtalResponse
      test case.
  */
  class TestCfg {
  public:
    TestCfg() {
      clear();
    }

    /// set all values to default
    void clear() {
      testXtals.clear();
      xtalHits.clear();
      trigMode      = idents::CalXtalId::BESTRANGE;
      zeroSuppress  = true;
    }

    /// select random test cfg
    void randomFill();

    /// max number of xtals to fill in one event
    static const unsigned short MAX_N_XTALS = 10;

    
    /// currently tested xtals
    typedef std::vector<CalUtil::XtalIdx> XtalList;

    /// currently tested xtals
    XtalList testXtals;

    /// count number of hit xtals in installed towers
    unsigned short nValidXtals(const TwrSet &twrSet) const;

    /// count installed xtals over LAC threshold plus some tolerance
    /// return < 0 on error
    short countXtalsOverLAC(ICalCalibSvc &calCalibSvc,
                            const TwrSet &twrSet,
                            float toleranceMeV) const;
  


    /// max mev in single xtal (just above HEX1 saturation
    static const float MAX_MEV;

    /// represent one single xtal hit
    struct XtalHit {
      XtalHit(const float _meV,
              const float _xtalPos) :
        meV(_meV),
        xtalPos(_xtalPos)
      {}

      /// energy of hit
      float meV;

      /// current test deposit position in fraction of crystal
      /// length. (from 0.0 -> 1.0)
      float xtalPos;
    };

    /// max number of mc hits for one event. (should never be more
    /// than number of MC segments per xtal
    static const unsigned short MAX_HITS_PER_XTAL = 2;

    /// collect all hits in one event.
    typedef std::vector<XtalHit> HitList;

    /// collect all hits in one event.
    HitList xtalHits;

    /// sum meV of all hits in single crystal
    float totalEnergy() const;

    /// weighted average energy centroid of hit
    float weightedPos() const;
    
    /// cal trig mode (BESTRANGE / ALLRANGE) for current test
    idents::CalXtalId::CalTrigMode trigMode;

    /// zero suppression for current test
    bool zeroSuppress;

  };

  ostream &operator<<(ostream &stream, const TestCfg &testCfg);
  
  ostream &operator<<(ostream &stream, const TestCfg::XtalHit &xtalHit);

}
#endif
