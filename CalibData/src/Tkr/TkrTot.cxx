// $Header$

/**
  @file TkrTot.cxx

  Implementation for tracker charge injection ToT, organized by
  (uniplane) Si layer.

  Implementation for two classes may be found here.  @a TkrTot is a class
  for per uniplane information.  @a TkrTotCol contains TkrTot
  information for the entire detector.
*/
#include "CalibData/Tkr/TkrTot.h"

namespace CalibData {

  /// Define factory class to TkrTotUni
  class UniTotFactory : public UniFactoryBase {
  public:
    UniTotFactory() : UniFactoryBase() {}
    virtual ~UniTotFactory() {}
    virtual UniBase* makeUni() {return new TkrTotUni();}
  };

  // Remainder of trivial TkrTotStrip class implementation
  void TkrTotStrip::copy(const TkrTotStrip& other) {
    m_stripId = other.m_stripId;
    m_slope = other.m_slope;
    m_intercept = other.m_intercept;
    m_quad = other.m_quad;
    m_chi2 = other.m_chi2;
    m_df = other.m_df;
  }
  void TkrTotStrip::init(int id, float slope, float intercept,
                         float quad, float chi2, float df) {
    m_stripId = id;
    m_slope = slope;
    m_intercept = intercept;
    m_quad = quad;
    m_chi2 = chi2;
    m_df = df;
  }


  // TkrTotUni
  TkrTotUni::TkrTotUni(const idents::TkrId& id, int nStrips) : 
    UniBase(id), m_nStrips(nStrips) {
    m_strips = new TkrTotStrip[nStrips];
  }

  void TkrTotUni::update(UniBase* other) {
    TkrTotUni* otherUni = dynamic_cast<TkrTotUni* > (other);
    using idents::TkrId;

    if (otherUni ==  0) return;

    if (otherUni->m_nStrips != m_nStrips) {
      m_nStrips = otherUni->m_nStrips;
      if (m_strips != 0) delete m_strips;
      m_strips = new TkrTotStrip[m_nStrips];
    }
    m_id = TkrId(otherUni->m_id);
    for (unsigned i = 0; i < m_nStrips; i++) {
      m_strips[i].copy(otherUni->m_strips[i]);
    }
  }

  bool TkrTotUni::putStrip(const TkrTotStrip& strip) {
    int iStrip = strip.getStripId();
    if ( (iStrip < 0) || (iStrip >= m_nStrips) ) return false;

    m_strips[iStrip].copy(strip);
    return true;
  }
  //  TkrTotCol

  TkrTotCol::TkrTotCol(unsigned nTowerRow, unsigned nTowerCol, unsigned nTray) 
    :    TkrBase(nTowerRow, nTowerCol, nTray) { 
    if (m_factory) delete m_factory;  // don't want base class factory
    m_factory = new UniTotFactory();
  }

  const CLID&  TkrTotCol::classID() {return CLID_Calib_TKR_TOTSignal;}

  bool TkrTotCol::putUni(UniBase* data, const idents::TkrId& id) {
    if (!dynamic_cast<TkrTotUni*>(data)) return false;

    // Otherwise go ahead and let base class handle it. 
    // We (well, TkrTotUni) might get called back to allocate memory
    return TkrBase::putUni(data, id);
  }

  /* maybe we don't need this after all
  StatusCode TkrTotCol::update(TkrBase& other, MsgStream* log) {

  }
  */


}
