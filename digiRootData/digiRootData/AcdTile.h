
#ifndef AcdTile_H
#define AcdTile_H

#include "AcdId.h"
#include "TObject.h"

/*! \class AcdTile
\brief 
 The digitization for a single ACD tile.
 Stores the pulse height and a bit to denote whether or not
 the tile is above veto threshold.
 Jun 2001 Heather Kelly renamed to AcdTile
 Jun 2001 Heather Kelly - Remove pointer to TileID
 Feb 2000 Daniel Flath - Rewrite to facilitate standardization of
			 code & interface across all classes.
 Jan 2000 Daniel Flath - ROOT HTML comments added
 Dec 1999 Daniel Flath - Rewritten for GLAST
 Oct 25,1999 Richard Dubois Clone from LCD version
*/

class AcdTile: public TObject {
private:
    /*!
    AcdTile tag word layout
    (defined in following enum):
     _______________________________________________
    |15|  |  |12|11|  |  |  |  |  |  |	|  |  |  |00|
    |__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
    |Hit Status |	      PHA Value 	    |
    |___________|___________________________________| 

    */

    enum {
	ACD_K_PMT = 12,
	ACD_V_PMT = 0,
	ACD_M_PMT = ((1 << ACD_K_PMT) - 1),
	ACD_K_HIT = 4,
	ACD_V_HIT = (ACD_K_PMT + ACD_V_PMT),
	ACD_M_HIT = ((1 << ACD_K_HIT) - 1)
    };
    /// packed word containing ACD tile data
    UShort_t m_tag;	
    /// Tile Id
    AcdId m_tileId; 
public:
    AcdTile();
    AcdTile(UInt_t id, short base=10, short used=1);
    AcdTile(AcdId &id);
    virtual ~AcdTile();
    AcdId* getId() { return &m_tileId; };
    UShort_t getPulseHeight();
    UChar_t getVeto();
    Bool_t setPulseHeight(UShort_t phaVal);
    Bool_t setVeto(UChar_t hitVal);

    /// Root >= 3.0 is now const correct for the Compare function
    Int_t Compare(const TObject *obj) const; 
    Bool_t IsSortable() const;

    ClassDef(AcdTile,3) 	// Digitization for a single ACD Tile
};

#endif
