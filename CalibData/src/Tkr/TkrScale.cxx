// $Header$

/**
  @file TkrScale.cxx

  Implementation for tracker charge scale, organized by
  (uniplane) Si layer.

  Implementation for three classes may be found here.  @a TkrScaleObj
  is a class is a trivial data class which has calib. information for one
  strip or gtfe.   @a TkrScale is a class
  for per uniplane information.  @a TkrScaleCol contains TkrScale
  information for the entire detector.
*/
#include "CalibData/Tkr/TkrScale.h"
#define N_STRIPS_PER_GTFE  64

namespace CalibData {

  /// Define factory class to TkrScaleUni
  class UniScaleFactory : public UniFactoryBase {
  public:
    UniScaleFactory() : UniFactoryBase() {}
    virtual ~UniScaleFactory() {}
    virtual UniBase* makeUni() {return new TkrScaleUni();}
  };

  TkrScaleCol::TkrScaleCol(unsigned nTowRow, unsigned nTowCol, unsigned nTray) 
    :    TkrBase(nTowRow, nTowCol, nTray) { 
    if (m_factory) delete m_factory;  // don't want base class factory
    m_factory = new UniScaleFactory();  
  }

  // Remainder of trivial TkrScaleObj class implementation
  void TkrScaleObj::copy(const TkrScaleObj& other) {
    m_id = other.m_id;
    m_scale = other.m_scale;
    m_error = other.m_error;
    m_chi2 = other.m_chi2;
    m_df = other.m_df;
  }
  void TkrScaleObj::init(int id, float scale, float error,
                         float chi2, float df) {
    m_id = id;
    m_scale = scale;
    m_error = error;
    m_chi2 = chi2;
    m_df = df;
  }


  // TkrScaleUni


  TkrScaleUni::TkrScaleUni(const idents::TkrId& id, int nObjs, 
                           bool perStrip) :
    UniBase(id), m_nObjs(nObjs), m_perStrip(perStrip) {
    m_objs = new TkrScaleObj[nObjs];
    m_nStrips = (perStrip) ? nObjs : nObjs * N_STRIPS_PER_GTFE;
  }

  void TkrScaleUni::update(UniBase* other) {
    TkrScaleUni* otherUni = dynamic_cast<TkrScaleUni* > (other);
    using idents::TkrId;

    if (otherUni ==  0) return;

    if (otherUni->m_nObjs != m_nObjs) {
      m_nObjs = otherUni->m_nObjs;
      if (m_objs != 0) delete m_objs;
      m_objs = new TkrScaleObj[m_nObjs];
    }
    m_perStrip = otherUni->m_perStrip;

    m_id = TkrId(otherUni->m_id);
    for (unsigned i = 0; i < m_nObjs; i++) {
      m_objs[i].copy(otherUni->m_objs[i]);
    }
  }

  bool TkrScaleUni::putObj(const TkrScaleObj& obj) {
    int id = obj.getId();
    if ( (id < 0) || (id >= m_nObjs) ) return false;

    m_objs[id].copy(obj);
    return true;
  }

  void TkrScaleUni::clearObjs() {
    for (unsigned ix=0; ix < m_nObjs; ix++) {
      (m_objs + ix)->clear();
    }
  }

  void TkrScaleUni::resize(unsigned n) {
    if (n == m_nObjs) clearObjs();
    else {
      delete [] m_objs;
      m_objs = new TkrScaleObj[n];
      m_nObjs = n;
    }
  }

  //  TkrScaleCol

  const CLID&  TkrScaleCol::classID() {return CLID_Calib_TKR_ChargeScale;}

  bool TkrScaleCol::putUni(UniBase* data, const idents::TkrId& id) {
    if (!dynamic_cast<TkrScaleUni*>(data)) return false;

    // Otherwise go ahead and let base class handle it. 
    // We (well, TkrScaleUni) might get called back to allocate memory
    return TkrBase::putUni(data, id);
  }
}
