//////////////////////////////////////////////////////////////////////
//
// The StripID class contains strip identification information.  This
// includes the strip number within the layer, the tower (always 0 for
// the '99/2000 beamtest, and the controller the strip was read out on.
//
//////////////////////////////////////////////////////////////////////

#include "digiRootData/StripID.h"

ClassImp(StripID)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
StripID::StripID() : m_tag(0)
{

}

/////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
Int_t StripID::Compare(TObject *obj) {
  StripID *strip;
  if (strip = dynamic_cast<StripID*>(obj)) {
    if (strip == this)
      return 0;
    else {
      UInt_t id_this = getID();
      UInt_t id_strip = strip->getID();
      if (id_this == id_strip)
	return 0;
      else
	return (id_this > id_strip) ? 1 : -1;
    }
  }
  else {
    return 0;
  }
}
/////////////////////////////////////////////////
Bool_t StripID::IsSortable() const {
  return kTRUE;
}
/////////////////////////////////////////////////
Bool_t StripID::setStrip(UShort_t stripVal) {
    if (!(stripVal & ~TKR_M_STRIP)) {
        m_tag &= ~(TKR_M_STRIP << TKR_V_STRIP);
        m_tag |= (stripVal << TKR_V_STRIP);
        return kTRUE;
    }
    else
        return kFALSE;
}

/////////////////////////////////////////////////
Bool_t StripID::setCTRL(UShort_t ctrlVal) {
    if (!(ctrlVal & ~TKR_M_CTRL)) {
        m_tag &= ~(TKR_M_CTRL << TKR_V_CTRL);
        m_tag |= (ctrlVal << TKR_V_CTRL);
        return kTRUE;
    }
    else
        return kFALSE;
}

/////////////////////////////////////////////////
Bool_t StripID::setTower(UShort_t towerVal) {
    if (!(towerVal & ~TKR_M_TOWER)) {
        m_tag &= ~(TKR_M_TOWER << TKR_V_TOWER);
        m_tag |= (towerVal << TKR_V_TOWER);
        return kTRUE;
    }
    else
        return kFALSE;
}
