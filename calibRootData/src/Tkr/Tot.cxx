
#include "calibRootData/Tkr/Tot.h"

ClassImp(calibRootData::TotStrip)

namespace calibRootData {

  // Return to state resulting from default constructor 
  void TotStrip::Clear(Option_t *) {
    m_stripId = -1;
    m_slope = 0.0;
    m_intercept = 0.0;
    m_quad = 0.0;
    m_chi2 = 0.0;
    m_df = 0.0;
  }

  void TotStrip::init(UInt_t id, Float_t slope, Float_t intercept,
                      Float_t quad, Float_t chi2, Float_t df) {
    m_stripId = id;
    m_slope = slope;
    m_intercept = intercept;
    m_quad = quad;
    m_chi2 = chi2;
    m_df = df;
  }

  void TotStrip::copy(const TotStrip& other) {
    m_stripId = other.m_stripId;
    m_slope = other.m_slope;
    m_intercept = other.m_intercept;
    m_quad = other.m_quad;
    m_chi2 = other.m_chi2;
    m_df = other.m_df;
  }
}

ClassImp(calibRootData::TotUnilayer)

namespace calibRootData {
  using commonRootData::TkrId;
  TotUnilayer::TotUnilayer(const TkrId& id, UInt_t nStrips) : 
    m_nStrips(nStrips), m_strips(0) {
    m_unilayerId = id;
    if (nStrips > 0) {
      m_strips = new TClonesArray("calibRootData::TotStrip", nStrips);
    }
  }

  void TotUnilayer::resize(unsigned nStrips) {
    if (nStrips == m_nStrips) return;
    //    if (m_strips) delete m_strips;
    if (m_strips) m_strips->Delete();
    m_nStrips = 0;
    if (nStrips == 0) return;
    m_strips = new TClonesArray("calibRootData::TotStrip", nStrips);
    if (m_strips) {
      m_nStrips = nStrips;
    }
  }

  Bool_t TotUnilayer::setStrip(const TotStrip& strip) {
    Int_t id = strip.getStripId();
    Int_t lastIx = m_strips->LastIndex();
    if (id > lastIx ) return false;
    //    TotStrip* pStrip = (TotStrip*) m_strips->At(id);
    //    pStrip->copy(strip);
    //    if (m_strips->At(id) != 0) {  // delete old
      //      delete (*m_strips)[id];
    //    }
    // extra argument to new is location where pointer to new memory
    // is to be stored
    new((*m_strips)[id]) TotStrip(strip);
    //    *((*m_strips)[id]) = strip;
    return true;
  }

  // Leave size of TCloneArray the same; otherwise clear out everything
  void TotUnilayer::Clear(Option_t *) {
    m_unilayerId = TkrId();

    m_strips->Clear("C");
    /*
    Int_t lastIx = m_strips->GetLast();
    if (lastIx < 0) return;
    UInt_t lastIxU = lastIx;
    for (UInt_t ix = 0; ix <= lastIxU; ix++) {

      TObject* pObj = m_strips->At(ix);
      TotStrip* pStrip = (TotStrip* ) pObj;
      if (pStrip) pStrip->Clear();
    }
    */
  }
}
