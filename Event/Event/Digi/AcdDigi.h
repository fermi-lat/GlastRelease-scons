#ifndef Event_AcdDigi_H
#define Event_AcdDigi_H 1

#include <iostream>
#include "idents/AcdId.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"

#include "Event/TopLevel/Definitions.h"

/** @class AcdDigi        
* @brief AcdDigi represents the output from one ACD entity.  An entity could
* be an ACD tile or a fiber.  Each of which would contain 2 PMTs.
* There are PHA and discriminator values for each PMT.
* Thus, there each member variable is an array of two entries.
* - Low Discriminator enables the PHA value
* - Veto Discriminator nominal ACD veto signal
* - High Discriminator is for CAL calibration - CNO
*             
* There are no set methods in this class, users are expected to fill
* the data members through the constructor.
* $Header$
*/

extern const CLID& CLID_AcdDigi;

namespace Event {
    class AcdDigi : virtual public ContainedObject  { 
        
    public:
        
        typedef enum {
            A = 0,
            B = 1
        } PmtId;
                
        AcdDigi(const idents::AcdId &id, unsigned short *pha, bool *veto, 
            bool *lowThresh, bool *highThresh) : m_id(id) 
        {  
            m_pulseHeight[0] = pha[0]; m_pulseHeight[1] = pha[1];
            m_veto[0] = veto[0]; m_veto[1] = veto[1];
            m_low[0] = lowThresh[0]; m_low[1] = lowThresh[1];
            m_high[0] = highThresh[0]; m_high[1] = highThresh[1];
        };
        
        virtual ~AcdDigi() { };
        
        /// Retrieve reference to class definition structure
        virtual const CLID& clID() const   { return AcdDigi::classID(); }
        static const CLID& classID()       { return CLID_AcdDigi; }
        
        /// Retrieve ACD identifier
        inline const idents::AcdId getId() const { return m_id; };
        
        /// Retrieve pulse height from one PMT
        inline unsigned short getPulseHeight(PmtId id) const { return m_pulseHeight[id]; };
        
        inline bool getVeto(PmtId id) const { return m_veto[id]; };
        
        inline bool getLowDiscrim(PmtId id) const { return m_low[id]; };
        
        inline bool getHighDiscrim(PmtId id) const { return m_high[id]; };
        
        
        /// Serialize the object for writing
        virtual StreamBuffer& serialize( StreamBuffer& s ) const;
        /// Serialize the object for reading
        virtual StreamBuffer& serialize( StreamBuffer& s );
        
        friend std::ostream& operator << ( std::ostream& s, const AcdDigi& obj ) {
            return obj.fillStream(s);
        }
        
        /// Fill the ASCII output stream
        virtual std::ostream& fillStream( std::ostream& s ) const;
        
        
    private:
        
        /// Acd ID
        idents::AcdId        m_id;
        /// pulse height
        unsigned short       m_pulseHeight[2];
        /// nominal Acd veto signal
        bool                 m_veto[2];
        /// 1 bit Low threshold discriminator - enables the PHA
        bool                 m_low[2];
        /// 1 bit High threshold discriminator - used for calibration of the CAL
        bool                 m_high[2];
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
            << m_id << ", "
            << "\n     PMT A pulse height    = "
            << EventFloatFormat( EventFormat::width, EventFormat::precision )
            << m_pulseHeight[0] << ", "
            << "      PMT A Discriminators (low, veto, high)   = ( "
            << m_low[0]   << " , "
            << m_veto[0] << " , "
            << m_high[0] << " )"
            << "\n     PMT B pulse height  = "
            << EventFloatFormat( EventFormat::width, EventFormat::precision )
            << m_pulseHeight[1] << ", "
            << "      PMT B Discriminators (low, veto, high)   = ( "
            << m_low[1]   << " , "
            << m_veto[1] << " , "
            << m_high[1] << " )\n";
        
    }
    
    
}// namespace Event

#endif    // Event_AcdDigi_H

