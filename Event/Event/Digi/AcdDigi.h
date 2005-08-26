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

/** @class AcdDigi        
* @brief AcdDigi represents the digitization output from one ACD entity.  
*
* An entity could be an ACD tile or a fiber.  Each of which would contain 
* 2 PMTs.  There are PHA and discriminator values for each PMT.
* Thus, there each member variable is an array of two entries.
* - Low Discriminator enables the PHA value
* - Veto Discriminator nominal ACD veto signal
* - High Discriminator is for CAL calibration - CNO
* So the AcdDigi is comprised of:
* - AcdId
* - Energy in MeV - as a check on the pha values
* - 2 Pulse Height values
* - 2 Low discriminators
* - 2 Veto discriminators
* - 2 High discriminators
*             
* There are no set methods in this class, users are expected to fill
* the data members through the constructor.
*
* @author Heather Kelly
* $Header$
*/

static const CLID& CLID_AcdDigi = InterfaceID("AcdDigi", 1, 0);

namespace Event {
    class AcdDigi : virtual public ContainedObject  { 
        
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

        // This keeps Gaudi happy
        AcdDigi() {}
                
        AcdDigi(const idents::AcdId &id, const idents::VolumeIdentifier &volId,
            double energy, unsigned short *pha, 
            bool *veto, bool *lowThresh, bool *highThresh) 
            : m_id(id), m_tileName(""), m_tileNumber(-1), m_volId(volId), m_energy(energy)
        {  
            m_pulseHeight[0] = pha[0]; m_pulseHeight[1] = pha[1];
            m_veto[0] = veto[0]; m_veto[1] = veto[1];
            m_low[0] = lowThresh[0]; m_low[1] = lowThresh[1];
            m_high[0] = highThresh[0]; m_high[1] = highThresh[1];
            m_range[0] = LOW; m_range[1] = LOW;
            m_error[0] = NOERROR; m_error[1] = NOERROR;
        };
        
        virtual ~AcdDigi() { };

        void initLdfParameters(const char* tileName, int tileNumber, Range *rangeVals, ParityError *errorVals) {
            m_tileName = tileName;
            m_tileNumber = tileNumber;
            m_range[0] = rangeVals[0]; m_range[1] = rangeVals[1];
            m_error[0] = errorVals[0]; m_error[1] = errorVals[1];
            m_error[2] = errorVals[2]; m_error[3] = errorVals[3];
        }
        
        /// Retrieve reference to class definition structure
        virtual const CLID& clID() const   { return AcdDigi::classID(); }
        static const CLID& classID()       { return CLID_AcdDigi; }
        
        /// Retrieve ACD identifier
        inline const idents::AcdId getId() const { return m_id; };
        inline const idents::VolumeIdentifier getVolId() const { return m_volId; };

        inline const char* getTileName() const { return m_tileName; };

        inline int getTileNumber() const { return m_tileNumber; };
        
        inline double getEnergy() const { return m_energy; };

        /// Retrieve pulse height from one PMT
        inline unsigned short getPulseHeight(PmtId id) const { return m_pulseHeight[id]; };
        
        inline bool getVeto(PmtId id) const { return m_veto[id]; };
        inline bool getHitMapBit(PmtId id) const { return m_veto[id]; };
        
        inline bool getLowDiscrim(PmtId id) const { return m_low[id]; };
        inline bool getAcceptMapBit(PmtId id) const { return m_low[id]; };
        
        inline bool getHighDiscrim(PmtId id) const { return m_high[id]; };
        inline bool getCno(PmtId id) const { return m_high[id]; };
        
        inline Range getRange(PmtId id) const { return m_range[id]; };

        inline ParityError getParityError(PmtId id) const { return m_error[id]; };
        inline ParityError getOddParityError(PmtId id) const { return m_error[id]; };
        inline ParityError getHeaderParityError(PmtId id) const { return m_error[id+2]; };

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
        const char*                m_tileName;
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
        /// referred to as the zero suppression bits
        bool                 m_low[2];
        /// 1 bit High threshold discriminator - used for calibration of the CAL
        /// used for calibration of the CAL
        /// This really should not be stored here any longer, we only set CNO
        /// per FREE board, so it is not easy to know which PMT caused it.
        bool                 m_high[2];
        /// Range setting either LOW (0) or HIGH (1)
        Range                m_range[2];
        /// Stores the parity error bit from LDF:  NOERROR (0), ERROR (1)
        /// 0,1 corespond to odd parity for each PMT and 2,3 is the header
        /// parity found in the AEM header which coresponds to CMT/Data Error
        /// in the ACD ICD.  Each FREE board should have its own header parity bit
        ParityError          m_error[4];
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
            << m_id << ", " << m_tileName << " "
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
            << "     PMT B Range: " << m_range[1]
            << "\n   Odd Parity:  " << m_error[1]
            << "\n   Header Parity: " << m_error[3]; 
    }
    
    
}// namespace Event

#endif    // Event_AcdDigi_H

