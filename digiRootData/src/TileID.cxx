///////////////////////////////////////////////////////////////////////////
//                                                                        
// The TileID class contains a tile ID number.  
// TileID is used in the ACDTile class.
//                                                                        
///////////////////////////////////////////////////////////////////////////

// Apr 2000 Daniel Flath - Minor changes to function names, etc.
// Jan 1999 Daniel Flath - ROOT HTML Documentation added
// Dec 1999 Daniel Flath - Rewrite for GLAST
// Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID

#include "digiRootData/TileID.h"
#include "TClass.h"

ClassImp(TileID)

///________________________________________________________________________
TileID::TileID() {
  m_ID = 0;
}
///________________________________________________________________________
TileID::TileID(UShort_t id = 0) : m_ID(id) {
  // Creates a TileID object with a tile ID number of id.
  // (ID defaults to a value of 0 {zero} if no value is supplied.)
}
///________________________________________________________________________
UShort_t TileID::getID() const { 
  // Returns the tile ID number

  return m_ID;
}
///________________________________________________________________________
void TileID::setID(UShort_t newVal) { 
  // Sets the tile ID number to newVal

  m_ID = newVal;
}
///________________________________________________________________________
/// Implement the streamer ourselves for now...to take advantage
/// of schema evolution - and to allow our TBEvent class to handle
/// both old (<= Root v2.25) Root files, and new (>= Root v3.00) files
void TileID::Streamer(TBuffer &R__b)
{
   // Stream an object of class TileID.

   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
      if (R__v > 1) 
      { 
          TileID::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          return;
      }
      /// Old Versions
      TObject::Streamer(R__b);
      R__b >> m_ID;
   } else {
       TileID::Class()->WriteBuffer(R__b, this);
   }
}