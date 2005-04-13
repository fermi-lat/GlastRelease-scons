
#ifndef ROOT_AcdDigi_H
#define ROOT_AcdDigi_H

#include "commonRootData/idents/AcdId.h"
#include "TObject.h"
#include "TString.h"
#include "commonRootData/idents/VolumeIdentifier.h"

/** @class AcdDigi
 * @brief The digitization for a single ACD entity.  
 * 
 * This digitization could represent a 
 * tile or a ribbon.  Each ACD entity has 2 PMTs associated with 
 * it.  Each PMT produces a pulse height and 3 discriminator
 * bits:  low, veto, high.  The PMT values are stored in a packed
 * word.
 * 
 * $Header$
*/

class AcdDigi: public TObject {

private:
    /**
    AcdDigi packed word layout
    (defined in following enum):
    
    @verbatim
     ____________________________________________________
    |15|14  |13  |12 |11|  |  |  |  |  |  |  |  |  |  |  |
    |__|____|____|___|__|__|__|__|__|__|__|__|__|__|__|__|
    |  |High|Veto|Low|	      PHA Value 	         |
    |__|____|____|___|___________________________________| 
    @endverbatim


    ldf packed word contains range and parity error information
    @verbatim
     ____________________________________________________
    | | | | | |  |  |  |  |  |  |  |  |  | 1   | 0   |
    |_|_|_|_|_|__|__|__|__|__|__|__|__|__|_____|_____|
    | | | | | |  |  |  |  |  |  |  |  |  |Range|Error|
    |_|_| |_|________________________________________| 
    @endverbatim

    */

    enum {
	PMT_WIDTH = 12,  // width of field
	PMT_SHIFT = 0,   // shift from right
	PMT_MAX = ((1 << PMT_WIDTH) - 1), // word representing bits in this field
        LOW_WIDTH = 1,
	LOW_SHIFT = (PMT_WIDTH + PMT_SHIFT),
        VETO_WIDTH = 1,
	VETO_SHIFT = (LOW_WIDTH + LOW_SHIFT),
        HIGH_WIDTH = 1,
        HIGH_SHIFT = (VETO_WIDTH + VETO_SHIFT),
    };
  
    enum {
        ERROR_WIDTH = 1,
        ERROR_SHIFT = 0,
        RANGE_WIDTH = 1,
        RANGE_SHIFT = (ERROR_WIDTH + ERROR_SHIFT)
    };

public:
	/** @defgroup AcdDigiGroup AcdDigi End-User Interface */
	/*@{*/
    typedef enum {
        A = 0,
        B = 1
    } PmtId;

    typedef enum {
        LOW = 0,
        HIGH = 1
    } Range;

    typedef enum {
        NOERROR = 0,
        ERROR = 1
    } ParityError;
	/*@}*/

    
    AcdDigi();

    AcdDigi(const commonRootData::AcdId& id);
    
    AcdDigi(const commonRootData::AcdId& id, const commonRootData::VolumeIdentifier& volId, Float_t energy, 
        UShort_t *pha, Bool_t *veto,
        Bool_t *low, Bool_t *high);
    
    virtual ~AcdDigi() { };

    void initLdfParameters(const char *tileName, int tileNumber, Range *rangeCol, ParityError *errCol);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

	/** @ingroup AcdDigiGroup */
	/*@{*/
	/// Returns the Monte Carlo energy deposited in this ACD detector
    Float_t getEnergy() { return m_energy; };

    const commonRootData::AcdId& getId() const { return m_id; };
    const commonRootData::VolumeIdentifier& getVolId() const { return m_volId; };

    /// Returns the PHA value for the PMT requested
    UShort_t getPulseHeight (AcdDigi::PmtId pmt) const;
    /// Returns True/False denoting whether this ACD Tile's veto threshold bit is on.
    Bool_t getVeto(AcdDigi::PmtId pmt) const;
    /// Returns true/false denoting if this Acd PMT high discrim is on/off
    Bool_t getLowDiscrim(AcdDigi::PmtId pmt) const;
    /// Returns true/false denoting if this ACD Tile's high discriminator bit is on (CNO).
    Bool_t getHighDiscrim(AcdDigi::PmtId pmt) const;

    Int_t getTileNumber() const { return m_tileNumber; };
    const char* getTileName() const { return m_tileName.Data(); };
    Range getRange(AcdDigi::PmtId pmt) const;
    ParityError getParityError(AcdDigi::PmtId pmt) const;
	/*@}*/

    /// Root >= 3.0 is now const correct for the Compare function
    Int_t Compare(const TObject *obj) const; 
    Bool_t IsSortable() const;

private:

    void initPackedWord(PmtId pmt, UShort_t pha, Bool_t veto, Bool_t low, Bool_t high);

    void initPackedLdfWord(PmtId pmt, Range range, ParityError error);

    /// Energy deposited in MeV - provided as a check
    Float_t m_energy;
    /// packed word containing ACD digi data
    UShort_t m_packed[2];	
    /// ACD Id
    commonRootData::AcdId m_id; 
    /// Volume id for geometry
    commonRootData::VolumeIdentifier m_volId;

    // New data members for version 3 
    UShort_t m_packedLdf[2];
    Int_t m_tileNumber;
    TString m_tileName;

    ClassDef(AcdDigi,3) // Digitization for a single ACD entity
};

#endif
