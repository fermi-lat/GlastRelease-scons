
#include "calibRootData/Tkr/ChargeScale.h"

ClassImp(calibRootData::ChargeScaleObj)

namespace calibRootData {

  // Return to state resulting from default constructor 
  void ChargeScaleObj::Clear(Option_t *) {
    m_Id = -1;
    m_scale = 0.0;
    m_error = 0.0;
    m_chi2 = 0.0;
    m_df = 0.0;
  }

  void ChargeScaleObj::init(UInt_t id, Float_t scale, Float_t error,
                            Float_t chi2, Float_t df) {
    m_Id = id;
    m_scale = scale;
    m_error = error;
    m_chi2 = chi2;
    m_df = df;
  }

  void ChargeScaleObj::copy(const ChargeScaleObj& other) {
    m_Id = other.m_Id;
    m_scale = other.m_scale;
    m_error = other.m_error;
    m_chi2 = other.m_chi2;
    m_df = other.m_df;
  }
}

ClassImp(calibRootData::ChargeScaleUnilayer)

namespace calibRootData {
  using commonRootData::TkrId;
  ChargeScaleUnilayer::ChargeScaleUnilayer(const TkrId& id) : 
    m_nChildren(0), m_children(0) {
    m_unilayerId = id;
    /*
    if (nStrips > 0) {
      m_strips = new TClonesArray("calibRootData::ChargeScaleStrip", nStrips);
    }
    */
  }

  void ChargeScaleUnilayer::resize(unsigned nObj, bool childIsStrip) {
    m_childIsStrip = childIsStrip;
    if (nObj == m_nChildren) return;
    if (m_children) m_children->Delete();
    m_nChildren = 0;
    if (nObj == 0) return;
    m_children = new TClonesArray("calibRootData::ChargeScaleObj", nObj);
    if (m_children) {
      m_nChildren = nObj;
    }
  }

  Bool_t ChargeScaleUnilayer::setChild(const ChargeScaleObj& child) {
    Int_t id = child.getId();
    Int_t lastIx = m_children->LastIndex();
    if (id > lastIx ) return false;
    new((*m_children)[id]) ChargeScaleObj(child);
    return true;
  }

  // Leave size of TCloneArray the same; otherwise clear out everything
  void ChargeScaleUnilayer::Clear(Option_t *) {
    m_unilayerId = TkrId();

    m_children->Clear("C");
  }
}
