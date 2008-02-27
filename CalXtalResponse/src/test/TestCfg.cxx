// $Header$
/** @file
    @author Z.Fewtrell

*/

// LOCAL
#include "TestCfg.h"
#include "test_util.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST

// EXTLIB
#include "CLHEP/Random/RandFlat.h"

// STD
#include <cassert>
#include <cmath>
#include <memory>
#include <list>
#include <set>

using namespace std;
using namespace CalUtil;

/// Utilities
namespace {
  /// shoot random number up to limit with exponential distribution
  float shootExp(float min, float max) {
    assert(max>min);
    assert(min > 0);
    
    return exp(CLHEP::RandFlat::shoot(log(min),log(max)));
  }

  /// shoot random unsigned integer < n with exponential distribution
  unsigned shootExpUInt(unsigned n) {
      /// add 1 to guarantee that log() works, and then subtract it later
    return static_cast<unsigned>(floor(exp(CLHEP::RandFlat::shoot(std::log(1.),std::log(n+1.))))-1);
  }

  static std::auto_ptr<CalXtalResponse::TestCfg> _testCfg;
  
  /// intended LAC threshold in mev
  const static float LAC_MEV_INTENT = 2.0;

}

namespace CalXtalResponse {
  const float TestCfg::MAX_MEV = 1e5; /// 100 gev

  /// fill each cfg parm w/ random value
  void TestCfg::randomFill() {
    clear();

    /// select number of xtals
    const unsigned short nXtals = CLHEP::RandFlat::shootInt(MAX_N_XTALS+1);
    
    /// create list of randomly selected xtals, (avoid duplicates using std set)
	std::set<XtalIdx> xtalSet;
	while (xtalSet.size() < nXtals) 
		xtalSet.insert(CalUtil::XtalIdx(CLHEP::RandFlat::shootInt(CalUtil::XtalIdx::N_VALS)));
    testXtals = vector<XtalIdx>(xtalSet.begin(), xtalSet.end());

    const float totalMeV = shootExp(0.1, MAX_MEV);

    /// i want exponential distribution with the most likely answer near MAX, not near zero - so I subtract.
    const unsigned short nHits = MAX_HITS_PER_XTAL - shootExpUInt(MAX_HITS_PER_XTAL+1);

    /// build nHits hits in different xtal segments with correct total
    /// energy.
    /// const list of all xtal segments
    static const unsigned short nSegments = 12;
    static const unsigned short fullXtalSegmentList[nSegments] = {0,1,2,3,4,5,6,7,8,9,10,11};
    /// list of unused segments for this crystal.
    typedef std::vector<unsigned short> SegmentList;
    SegmentList remainingSegments(fullXtalSegmentList,
                                  fullXtalSegmentList+nSegments);
    float remainingEne = totalMeV;
    for (unsigned short n = 0;
         n < nHits;
         n++) {
      float hitEne(0);
      if (n == nHits-1)  // use all remaining energy
        hitEne = remainingEne;
      else
        /// this hit should have some fraction of total energy (use at least .1% to avoid ridiculously small deposits.
        hitEne = CLHEP::RandFlat::shoot(0.,remainingEne);
      
      remainingEne -= hitEne;

      /// pick random segment from list of unused segments
      const unsigned short segmentIdx = CLHEP::RandFlat::shootInt(remainingSegments.size());
      unsigned short segment = remainingSegments[segmentIdx];
      remainingSegments.erase(remainingSegments.begin() + segmentIdx);

      /// make sure each hit is in different segment of xtal
      const float xtalSegmentMin = (float)segment/(float)nSegments;
      const float xtalSegmentMax = (float)(segment+1)/(float)nSegments;
      const float xtalPos = CLHEP::RandFlat::shoot(xtalSegmentMin, xtalSegmentMax);

      xtalHits.push_back(XtalHit(hitEne, xtalPos));
    }

    trigMode = CLHEP::RandFlat::shootBit() ? idents::CalXtalId::BESTRANGE : idents::CalXtalId::ALLRANGE;


    zeroSuppress = CLHEP::RandFlat::shootBit();

    createDiagnosticData = CLHEP::RandFlat::shootBit();
      
  }

  float TestCfg::totalEnergy() const {
    float retVal = 0;
    
    for (HitList::const_iterator it(xtalHits.begin());
         it != xtalHits.end();
         it++) 
      retVal += it->meV;

    return retVal;
  }

  float TestCfg::weightedPos() const {
    float posSum = 0;
    float eneSum = 0;

    for (HitList::const_iterator it(xtalHits.begin());
         it != xtalHits.end();
         it++)  {
      posSum += it->xtalPos*it->meV;
      eneSum += it->meV;
    }

    return posSum / eneSum;
  }

  ostream &operator<<(ostream &stream, const TestCfg &testCfg) {
    stream << "CalXtalResponse test case:" 
           << " N_XTALS=" << testCfg.testXtals.size()
           << " N_HITS=" << testCfg.xtalHits.size()
           << " TRIG_MODE " << testCfg.trigMode
           << " zeroSuppress " << testCfg.zeroSuppress
		   << " diag " << testCfg.createDiagnosticData
           << endl;

    for (TestCfg::XtalList::const_iterator xtalIt(testCfg.testXtals.begin());
         xtalIt != testCfg.testXtals.end();
         xtalIt++) 
      stream << " XTAL=" << xtalIt->toStr() << endl;

    for (TestCfg::HitList::const_iterator hitIt(testCfg.xtalHits.begin());
         hitIt != testCfg.xtalHits.end();
         hitIt++)
      stream << " HIT=[" << *hitIt << "]" <<endl;

    return stream;
  }

  ostream &operator<<(ostream &stream, const TestCfg::XtalHit &xtalHit) {
    stream << " MeV=" << xtalHit.meV 
           << " Pos=" << xtalHit.xtalPos;
    
    return stream;
  }

  unsigned short TestCfg::nValidXtals(const TwrSet &twrSet) const {
    unsigned short retVal(0);
    for (TestCfg::XtalList::const_iterator xtalIt(testXtals.begin());
         xtalIt != testXtals.end();
         xtalIt++) 
      if (twrSet.find(xtalIt->getTwr()) != twrSet.end())
        retVal++;
    return retVal;
  }

  /// count installed xtals over LAC threshold plus some tolerance
  short TestCfg::countXtalsOverLAC(ICalCalibSvc &calCalibSvc,
                                   const TwrSet &twrSet,
                                   const float toleranceMeV) const {
    const float ene = totalEnergy();
    const float relPos = weightedPos();

    unsigned short xtalCount(0);
    for (TestCfg::XtalList::const_iterator xtalIt(testXtals.begin());
         xtalIt != testXtals.end();
         xtalIt++) 
      if (twrSet.find(xtalIt->getTwr()) != twrSet.end()) {
        const float posMM = xtalRelPos2MM(relPos);

        float asymLrgDiode = 0;
        if (calCalibSvc.evalAsym(*xtalIt, ASYM_LL, posMM, asymLrgDiode).isFailure())
          return -1;
        
        const float asymRatio = exp(asymLrgDiode);

        // here we go again
        // ene = sqrt(pos*neg)
        // pos = e^2/neg
        // neg = e^2/pos
        //
        // asymRatio = pos/neg
        // pos = neg*asym
        // neg = pos/asym
        //
        // pos = asym*e^2/pos
        // pos^2 = asym*e^2
        // pos = e*sqrt(asym)
        //
        // neg = e^2/(neg*asym)
        // neg^2 = e^2/asym
        // neg = e/sqrt(asym)
        
        const float posMeV = ene*sqrt(asymRatio);
        const float negMeV = ene/sqrt(asymRatio);
    
        if (posMeV > LAC_MEV_INTENT + toleranceMeV  ||
            negMeV> LAC_MEV_INTENT + toleranceMeV  )
          xtalCount++;
      }

    return xtalCount;
  }

                                       

} // namespace CalXtalResponse
