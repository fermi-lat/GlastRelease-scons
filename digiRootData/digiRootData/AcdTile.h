
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
    |15|  |13|12|11|  |  |  |  |  |  |	|  |  |  |00|
    |__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
    |      CNOVe|	      PHA Value 	    |
    |___________|___________________________________| 

    */

    enum {
	ACD_K_PMT = 12,  // width of field
	ACD_V_PMT = 0,   // shift from right
	ACD_M_PMT = ((1 << ACD_K_PMT) - 1), // word representing bits in this field
	ACD_K_VETO = 1,
	ACD_V_VETO = (ACD_K_PMT + ACD_V_PMT),
	ACD_M_VETO = ((1 << ACD_K_VETO) - 1),
        ACD_K_HIGH = 1,
        ACD_V_HIGH = (ACD_K_VETO + ACD_V_VETO),
        ACD_M_HIGH = ((1 << ACD_K_HIGH) - 1)
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
    /// Returns the PHA value for this ACD Tile
    UShort_t getPulseHeight();
    /// Returns 1/0 denoting whether this ACD Tile's veto threshold bit is on.
    UChar_t getVeto();
    /// Returns 1/0 denoting if this ACD Tile's high discriminator bit is on (CNO).
    UChar_t getCNO();

    Bool_t setPulseHeight(UShort_t phaVal);
    Bool_t setVeto(UChar_t hitVal);
    Bool_t setCNO(UChar_t hitVal);

    /// Root >= 3.0 is now const correct for the Compare function
    Int_t Compare(const TObject *obj) const; 
    Bool_t IsSortable() const;

    ClassDef(AcdTile,4) 	// Digitization for a single ACD Tile
};

#endif
