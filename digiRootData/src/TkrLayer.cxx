////////////////////////////////////////////////////////////////////////
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
//#include "TClass.h"

ClassImp(TkrLayer)

TObjArray *TkrLayer::m_staticArray = 0;

// -------------------------------------------------------------------
TkrLayer::TkrLayer() {
  // Default constructor
    if (!m_staticArray) m_staticArray = new TObjArray();
    m_strips = m_staticArray;
}
// -------------------------------------------------------------------
/*
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
*/
// -------------------------------------------------------------------
TkrLayer::~TkrLayer() {
  // Destructor.  Removes all strips from 

  if (m_strips) {
    //int nEntries = m_strips->GetEntries();
   // for (int i=0; i<nEntries; i++)
   //   delete m_strips->At(i);
    m_strips->Delete();
    if (m_strips == m_staticArray) m_staticArray = 0;
    delete m_strips;
    m_strips = 0;
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

void TkrLayer::Clean(Option_t *option) {
    m_strips->Delete(option);
   // delete m_strips;
 //   m_strips = 0;
}


