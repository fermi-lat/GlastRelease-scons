//////////////////////////////////////////////////////////////////////
//
// The TkrLayer class provides access to a single layer of Si Strips.
// It contains a list of strips that were above thresh, and TOT and 
// ERRF values.
//
//////////////////////////////////////////////////////////////////////

// --------------------- Modification History ------------------------
// Apr 2000 Daniel Flath - Cleaned up to look like other files in package
// Jan 1999 Daniel Flath - ROOT HTML Documentation added
// Dec 1999 Daniel Flath - Creation

#include "digiRootData/TkrLayer.h"
#include "TClass.h"

ClassImp(TkrLayer)

// -------------------------------------------------------------------
TkrLayer::TkrLayer() {
  // Default constructor
}
// -------------------------------------------------------------------
TkrLayer::TkrLayer(TObjArray *strips = 0) {
  // Creates a new TkrLayer object with the specified TObjArray of Si
  // Si strips.  
  // (The array pointer defaults to 0 {NULL}, in which case a new
  // TObjArray is created.)
  
  if (strips)
    m_strips = strips;
  else
    m_strips = new TObjArray();
}
// -------------------------------------------------------------------
TkrLayer::~TkrLayer() {
  // Destructor.  Removes all strips from 

  if (m_strips) {
    int nEntries = m_strips->GetEntries();
    for (int i=0; i<nEntries; i++)
      delete m_strips->At(i);
    m_strips->Clear();
    delete m_strips;
  }
}
// -------------------------------------------------------------------
Int_t TkrLayer::Compare(const TObject *obj) const {
    if (this == obj) return 0;
    if (TkrLayer::Class() != obj->IsA()) return -1;
    UInt_t id_this = (getLayer() << 1) | (getXY() ? 0 : 1);
    UInt_t id_layer = (((TkrLayer*)obj)->getLayer() << 1) | (((TkrLayer*)obj)->getXY() ? 0 : 1);
    if (id_this == id_layer)
        return 0;
    else
	return (id_this > id_layer) ? 1 : -1;
}
// -------------------------------------------------------------------
Bool_t TkrLayer::IsSortable() const {
  return kTRUE;
}
// -------------------------------------------------------------------
void TkrLayer::setTOT(UInt_t right, UInt_t left) {
    m_TOT[0] = right;
    m_TOT[1] = left;
}
// -------------------------------------------------------------------
void TkrLayer::setErrf(UInt_t right, UInt_t left) {
    m_errf[0] = right;
    m_errf[1] = left;
}
// -------------------------------------------------------------------
Int_t TkrLayer::getTOT(UChar_t ctrlNum) {
    if ((ctrlNum == 0) || (ctrlNum ==1))
        return m_TOT[ctrlNum];
    else
        return -1;
}
// -------------------------------------------------------------------
Int_t TkrLayer::getErrf(UChar_t ctrlNum) {
    if ((ctrlNum == 0) || (ctrlNum ==1))
        return m_errf[ctrlNum];
    else
        return -1;
}
// -------------------------------------------------------------------
/// Implement the streamer ourselves for now...to take advantage
/// of schema evolution - and to allow our TBEvent class to handle
/// both old (<= Root v2.25) Root files, and new (>= Root v3.00) files
void TkrLayer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TkrLayer.
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
      if (R__v > 1) 
      { 
          TkrLayer::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          return;
      }
      /// Old Versions
      TObject::Streamer(R__b);
      R__b >> m_layer;
      R__b >> m_xy;
      R__b.ReadStaticArray(m_TOT);
      R__b.ReadStaticArray(m_errf);
      R__b >> m_strips;
   } else {
       TkrLayer::Class()->WriteBuffer(R__b, this);
   }
}




