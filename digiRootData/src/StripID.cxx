//////////////////////////////////////////////////////////////////////
//
// The StripID class contains strip identification information.  This
// includes the strip number within the layer, the tower (always 0 for
// the '99/2000 beamtest, and the controller the strip was read out on.
//
//////////////////////////////////////////////////////////////////////

#include "digiRootData/StripID.h"
#include "TClass.h"

ClassImp(StripID)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
StripID::StripID() : m_tag(0)
{

}

/////////////////////////////////////////////////
StripID::~StripID()
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

/////////////////////////////////////////////////
/// Implement the streamer ourselves for now...to take advantage
/// of schema evolution - and to allow our TBEvent class to handle
/// both old (<= Root v2.25) Root files, and new (>= Root v3.00) files
void StripID::Streamer(TBuffer &R__b)
{
   // Stream an object of class StripID.

   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
      if (R__v > 1) 
      { 
          StripID::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          return;
      }
      /// Old Versions      
      TObject::Streamer(R__b);
      R__b >> m_tag;
   } else {
       StripID::Class()->WriteBuffer(R__b, this);
   }
}