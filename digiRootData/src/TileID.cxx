///////////////////////////////////////////////////////////////////////////
//									  
// The TileID class contains a tile ID number.	
// TileID is used in the AcdTile class.
//									  
///////////////////////////////////////////////////////////////////////////

// Apr 2000 Daniel Flath - Minor changes to function names, etc.
// Jan 1999 Daniel Flath - ROOT HTML Documentation added
// Dec 1999 Daniel Flath - Rewrite for GLAST
// Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID

#include "digiRootData/TileID.h"
//#include "TClass.h"

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
