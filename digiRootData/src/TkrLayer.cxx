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

ClassImp(TkrLayer)

// -------------------------------------------------------------------
TkrLayer::TkrLayer() {
  // Default constructor
    m_strips = 0;
    m_numStrips = -1;
}
// -------------------------------------------------------------------

TkrLayer::TkrLayer(TObjArray *strips = 0) {
  // Creates a new TkrLayer object with the specified TObjArray of Si
  // Si strips.  
  // (The array pointer defaults to 0 {NULL}, in which case a new
  // TObjArray is created.)
  m_numStrips = -1;
  if (strips)
    m_strips = strips;
  else
    m_strips = new TObjArray();
}

// -------------------------------------------------------------------
TkrLayer::~TkrLayer() {
  // Destructor.  Removes all strips
  if (m_strips) {
    m_strips->Delete();
    delete m_strips;
    m_strips = 0;
  }
}
// -------------------------------------------------------------------
Int_t TkrLayer::Compare(const TObject *obj) const {
    // Ritz numbering convention:  Layer Number 0 is the bottom most layer - closest to the CAL
    if (this == obj) return 0;
    if (TkrLayer::Class() != obj->IsA()) return -1;

    // Even numbered planes have the layer measuring Y on the bottom - closer to CAL
    // Odd numbered planes have the layer measuring X on the bottom - closer to CAL
    UInt_t id_this = getLayerNum();
    UInt_t id_layer = ((TkrLayer*)obj)->getLayerNum();

    //id_this = (getPlaneNum() << 1) | (getXY() ? 0 : 1);    
    //id_layer = (((TkrLayer*)obj)->getPlaneNum() << 1) | (((TkrLayer*)obj)->getXY() ? 0 : 1);
    
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
void TkrLayer::setToT(UInt_t right, UInt_t left) {
    m_ToT[0] = right;
    m_ToT[1] = left;
}
// -------------------------------------------------------------------
void TkrLayer::setErrf(UInt_t right, UInt_t left) {
    m_errf[0] = right;
    m_errf[1] = left;
}
// -------------------------------------------------------------------
Int_t TkrLayer::getToT(UChar_t ctrlNum) {
    if ((ctrlNum == 0) || (ctrlNum ==1))
	return m_ToT[ctrlNum];
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

UShort_t TkrLayer::getLayerNum() const {
    // Even numbered planes have the layer measuring Y on the bottom - closer to CAL
    // Odd numbered planes have the layer measuring X on the bottom - closer to CAL
    if ( (getPlaneNum() % 2) == 0 ) {   // even plane
        // So in this case, we want the layer measuring Y to have the smaller id
        return ( (getPlaneNum() << 1) | (getXY() ? 0 : 1) );
    } else  {                           // odd plane
        // here we want the layer measuring X to have the smaller id
        return ( (getPlaneNum() << 1) | (getXY() ? 1 : 0) );
    }
}

void TkrLayer::Clean(Option_t *option) {
    m_strips->Delete(option);
    m_numStrips = -1;
   // delete m_strips;
 //   m_strips = 0;
}

