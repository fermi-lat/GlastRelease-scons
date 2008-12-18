#ifndef Event_AcdDigi_H
#define Event_AcdDigi_H 1

#include <iostream>
#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/IInterface.h"

#include "Event/TopLevel/Definitions.h"

/** @class Event::AcdDigi        
* @brief AcdDigi represents the digitization output from one ACD entity.  
*
* An entity could be an ACD tile or a fiber.  Each of which would contain 
* 2 PMTs.  There are PHA and discriminator values for each PMT.
* Thus, there each member variable is an array of two entries.
* - Low Discriminator (Accept Map bit == zero suppression) enables the PHA value
* - Veto Discriminator (Hit Map bit) nominal ACD veto signal
* - High (CNO) Discriminator is for CAL calibration
* So the AcdDigi is comprised of:
* - AcdId
* - Energy in MeV - as a check on the pha values
* - 2 Pulse Height values
* - 2 Low discriminators (Accept Map)
* - 2 Veto discriminators (Hit Map)
* - 2 High (CNO) discriminators, really only used for MC, real data has CNO
*   stored in GEM.  AcdDigi alg needs to be updated so MC behaves the same
*             
*
* @author Heather Kelly
* $Header$
*/

static const CLID& CLID_AcdDigi = InterfaceID("AcdDigi", 1, 1);

namespace Event {

class AcdDigi : virtual public ContainedObject  
{ 
public:
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

    /// Status bit definitions
    enum DigiStatusBits {
         DIGI_NULL    = 0,
         DIGI_MC      = 1<<28,   // This digi is from Monte Carlo
         DIGI_DATA    = 1<<30,   // This digi is from data 
         DIGI_OVERLAY = 1<<31    // This digi contains overlay info
    };

    // This keeps Gaudi happy
    AcdDigi() {}
            
    AcdDigi(const idents::AcdId &id, const idents::VolumeIdentifier &volId, const char *tileName="")
      : m_id(id), m_tileName(tileName), m_tileNumber(-1), m_volId(volId),
        m_energy(0.0f), m_ninja(false), m_gem(false), m_status(DIGI_NULL) 
    {
        m_pulseHeight[0] = 0; m_pulseHeight[1] = 0;
        m_veto[0] = false; m_veto[1] = false;
        m_low[0] = false; m_low[1] = false;
        m_high[0] = false; m_high[1] = false;
        m_range[0] = LOW; m_range[1] = LOW;
        m_error[0] = NOERROR; m_error[1] = NOERROR;
        m_error[2] = NOERROR; m_error[3] = NOERROR;
    }


    AcdDigi(const idents::AcdId &id, const idents::VolumeIdentifier &volId,
        double energy, unsigned short *pha, 
        bool *veto, bool *lowThresh, bool *highThresh) 
        : m_id(id), m_tileName(""), m_tileNumber(-1), m_volId(volId), 
          m_energy(energy), m_ninja(false), m_gem(false), m_status(DIGI_NULL)
    {  
        m_pulseHeight[0] = pha[0]; m_pulseHeight[1] = pha[1];
        m_veto[0] = veto[0]; m_veto[1] = veto[1];
        m_low[0] = lowThresh[0]; m_low[1] = lowThresh[1];
        m_high[0] = highThresh[0]; m_high[1] = highThresh[1];
        m_range[0] = LOW; m_range[1] = LOW;
        m_error[0] = NOERROR; m_error[1] = NOERROR;
    m_error[2] = NOERROR; m_error[3] = NOERROR;
    };
    
    virtual ~AcdDigi() { };

    void init(unsigned short* pha, bool *veto, bool* lowThresh, bool* highThresh, unsigned int status=DIGI_NULL) 
    {
        m_pulseHeight[0] = pha[0]; m_pulseHeight[1] = pha[1];
        m_veto[0] = veto[0]; m_veto[1] = veto[1];
        m_low[0] = lowThresh[0]; m_low[1] = lowThresh[1];
        m_high[0] = highThresh[0]; m_high[1] = highThresh[1];
        m_status = status;
    }

    void initLdfParameters(const char* tileName, int tileNumber, Range *rangeVals, ParityError *errorVals) {
        m_tileName = tileName;
        m_tileNumber = tileNumber;
        m_range[0] = rangeVals[0]; m_range[1] = rangeVals[1];
        m_error[0] = errorVals[0]; m_error[1] = errorVals[1];
        m_error[2] = errorVals[2]; m_error[3] = errorVals[3];
    }

    void initGem(bool ninja, bool gem) {
        setNinja(ninja);
        setGem(gem);
    }

    void setNinja(bool val=true) { m_ninja = val; }

    void setGem(bool val=true) { m_gem = val; }

    /// Allow status word to be modified
    inline void setStatus(unsigned int status)       {m_status  = status;}
    inline void addToStatus(DigiStatusBits bitToAdd) {m_status |= bitToAdd;}

    /// Set the ranges.  This is used by the AcdDigiAlg to hack in the correct 
    /// ranges in simulated data
    void setRanges(Range *rangeVals) {
        m_range[0] = rangeVals[0]; m_range[1] = rangeVals[1];
    }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return AcdDigi::classID(); }
    static const CLID& classID()       { return CLID_AcdDigi; }
    
    /// Retrieve ACD identifier
    inline const idents::AcdId getId() const { return m_id; };
    inline const idents::VolumeIdentifier getVolId() const { return m_volId; };

    /// Retrieve string name of the ACD detector as reported by LDF
    /// Corresponds to the ACD ids we are familiar with
    inline const char* getTileName() const { return m_tileName; };

    /// Retrieve tile number as reported by LDF, which may not correspond
    /// to accepted ACD ids
    inline int getTileNumber() const { return m_tileNumber; };
    
    /// Returns MC energy (MeV) deposited in the detector
    /// strictly here as a check on MC digitization algorithm
    inline double getEnergy() const { return m_energy; };

    /// Retrieve pulse height from one PMT
    inline unsigned short getPulseHeight(PmtId id) const { return m_pulseHeight[id]; };
    
    /// deprecated method name, see getHitMapBit
    inline bool getVeto(PmtId id) const { return m_veto[id]; };
    /// Denotes that the PMT was above hit (veto) threshold
    inline bool getHitMapBit(PmtId id) const { return m_veto[id]; };
    
    /// deprecated method name, see getAcceptMapBit
    inline bool getLowDiscrim(PmtId id) const { return m_low[id]; };
    /// Denotes that the PMT's PHA was read out 
    inline bool getAcceptMapBit(PmtId id) const { return m_low[id]; };
    
    /// Only useful for MC data for now
    inline bool getHighDiscrim(PmtId id) const { return m_high[id]; };
    /// Only useful for MC data for now
    inline bool getCno(PmtId id) const { return m_high[id]; };
    
    inline Range getRange(PmtId id) const { return m_range[id]; };

    //inline ParityError getParityError(PmtId id) const { return m_error[id]; };
    /// Error bit associated with PHA
    inline ParityError getOddParityError(PmtId id) const { return m_error[id]; };
    /// Error bit stored in AEM header
    inline ParityError getHeaderParityError(PmtId id) const { return m_error[id+2]; };

    /// true if this AcdDigi was solely created due to GEM TileList bit for this detector
    /// being set
    inline bool isNinja() const { return m_ninja; }

    /// true if GEM TileList bit for this detector is set
    inline bool getGemFlag() const { return m_gem; }

    /// Retrieve the status word
    inline unsigned int getStatus() const {return m_status;}

    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    
    friend std::ostream& operator << (std::ostream& s, const AcdDigi& obj)
    {
        return obj.fillStream(s);
    };
    
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;
    

private:
    
    /// Acd ID
    idents::AcdId        m_id;
    /// Tile name in char* form, matches idents::AcdId
    const char*          m_tileName;
    /// Tile Number as reported from LDF
    int                  m_tileNumber;
    /// Allow one to retrieve dimensions of this volume
    idents::VolumeIdentifier m_volId;
    /// Monte Carlo energy MeV - storing as a check on pulseHeight
    double               m_energy;
    /// pulse height PHA
    unsigned short       m_pulseHeight[2];
    /// nominal Acd veto signal
    bool                 m_veto[2];
    /// 1 bit Low threshold discriminator - enables the PHA
    /// Now more properly called the accept map bits and sometimes
    // referred to as the zero suppression bits
    bool                 m_low[2];
    /// 1 bit High threshold CNO discriminator 
    /// used for calibration of the CAL
    /// This really should not be stored here any longer.. we only set
    /// CNO per FREE board.  This data member remains for the MC.
    bool                 m_high[2];
    /// Range setting either LOW (0) or HIGH (1)
    Range                m_range[2];
    /// Stores the parity error bits from LDF:  NOERROR (0), ERROR (1)
    /// 0,1 correspond to odd parity for each PMT and 2,3 is the header
    /// parity found in the AEM header which corresponds to CMD/Data Error
    // in ACD ICD.  Each FREE board should have its own header parity bit
    ParityError          m_error[4];
    /// Flag denoting that this AcdDigi was generated solely due to the
    /// GEM bit in GemTileList for this detector
    bool                 m_ninja;
    /// Denotes that the bit corresponding to this detector in GemTileList
    /// is on.  Note that the GEM bits are detector level (not PMT specific)
    bool                 m_gem;
    /// Status bit word to keep track of source/status/etc.
    unsigned int         m_status;
};


//! Definition of all container types of AcdDigi
typedef ObjectVector<AcdDigi>     AcdDigiCol;

/// Serialize the object for writing
inline StreamBuffer& AcdDigi::serialize( StreamBuffer& s ) const
{
    ContainedObject::serialize(s);
    return s << m_id;
}


/// Serialize the object for reading
inline StreamBuffer& AcdDigi::serialize( StreamBuffer& s )
{
    ContainedObject::serialize(s);
    
    unsigned int id;
    
    s >> id;
    
    idents::AcdId persId(id);
    m_id = persId;
    
    return s;
}


/// Fill the ASCII output stream
inline std::ostream& AcdDigi::fillStream( std::ostream& s ) const
{
    return s
        << "    base class AcdDigi :"
        << "\n        id = ( "
        << EventFloatFormat( EventFormat::width, EventFormat::precision )
        << m_id << ", " << m_tileName << "   "
        << "\n     PMT A pulse height    = "
        << EventFloatFormat( EventFormat::width, EventFormat::precision )
        << m_pulseHeight[0] << ", "
        << "      PMT A Discriminators (accept, veto, cno)   = ( "
        << m_low[0]   << " , "
        << m_veto[0] << " , "
        << m_high[0] << " )"
        << "\n     PMT A Range: " << m_range[0]
        << "\n     Odd Parity: " << m_error[0]
        << "\n     Header Parity: " << m_error[2]
        << "\n     PMT B pulse height  = "
        << EventFloatFormat( EventFormat::width, EventFormat::precision )
        << m_pulseHeight[1] << ", "
        << "      PMT B Discriminators (accept, veto, cno)   = ( "
        << m_low[1]   << " , "
        << m_veto[1] << " , "
        << m_high[1] << " )\n"
        << "      PMT B Range: " << m_range[1]
        << "\n     Odd Parity: " << m_error[1]
        << "\n     Header Parity: " << m_error[3];
}


}// namespace Event

#endif    // Event_AcdDigi_H

