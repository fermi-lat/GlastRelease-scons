// TileID.h

// Apr 2000 Daniel Flath - Minor changes to function names, etc.
// Jan 1999 Daniel Flath - ROOT HTML Documentation added
// Dec 1999 Daniel Flath - Rewrite for GLAST
// Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID

#ifndef TILEID_H
#define TILEID_H

#include "TMath.h"
#include "TObject.h"

class TileID : public TObject {
 private:
  UShort_t m_ID;   // ACD Tile ID number
 public:
  TileID();
  TileID(UShort_t id);
  ~TileID() {};
  
  UShort_t getID() const;
  
  void setID(UShort_t newVal);
  
  ClassDef(TileID,2)      // ACD Tile identification information
};

#endif
