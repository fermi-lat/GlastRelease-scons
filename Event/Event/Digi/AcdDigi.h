#ifndef GlastEvent_AcdDigi_H
#define GlastEvent_AcdDigi_H 1


// Include files
#include <iostream>
#include "idents/AcdId.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"

#include "Event/TopLevel/Definitions.h"

/*!
//------------------------------------------------------------------------------
//
\class   AcdDigi        
  
\brief Digitizations for ACD                                
              
Initial Implementation by Richard Dubois
Additions by Heather Kelly

Note the usage of AcdId, this will be replaced by
the new AcdId class under development.

There are no set methods in this class, users are expected to fill
the data members through the constructor.

ACDDigi represents the output from one phototube.  The bits are organized
as follows:
     __________________________________________________________
    |15 | 14   | 13   | 12 |11|  |  |  |  |  |  |  |  |  |  |00|
    |__ |_____ |__ ___|____|__|__|__|__|__|__|__|__|__|__|__|__|
    |   |  HI  | LOW  |VETO|           PHA Value               |
    |___|THRESH|THRESH|____|___________________________________| 

  Low Threshold enables the PHA value
  Veto Threshold nominal ACD veto signal
  High Threshold is for CAL calibration - CNO

//------------------------------------------------------------------------------
 */

extern const CLID& CLID_AcdDigi;


class AcdDigi : virtual public ContainedObject  { 

public:
    AcdDigi() {};

    AcdDigi(idents::AcdId id, int pulseHeight, bool veto, bool lowThresh, 
        bool highThresh)       
        : m_ID(id),
        m_pulseHeight(pulseHeight),
        m_veto(veto),
        m_lowThreshold(lowThresh),
        m_highThreshold(highThresh),
        m_packedPHA(0)
    {  
        m_packedPHA = m_pulseHeight;  // set the lower 12 bits
        m_packedPHA |= (m_lowThreshold << ADC_V_LOWTHRESH);
        m_packedPHA |= (m_veto << ADC_V_VETO);
        m_packedPHA |= (m_highThreshold << ADC_V_HIGHTHRESH);
    };

    /// Destructor
    virtual ~AcdDigi() { };

    //! Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return AcdDigi::classID(); }
    static const CLID& classID()       { return CLID_AcdDigi; }

    /// Retrieve digi identifier
    inline const idents::AcdId ID() const { return m_ID; };

    /// Retrieve pulse height
    inline int PulseHeight() const { return m_pulseHeight; };

    inline bool Veto() const { return m_veto; };
    
    inline bool LowThreshold() const { return m_lowThreshold; };

    inline bool HighThreshold() const { return m_highThreshold; };


    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;



private:


    enum {
        ADC_K_VETO = 1,
        ADC_V_VETO = 12,    // determines location in packed byte
        ADC_M_VETO = ((1 << ADC_K_VETO) - 1),

        ADC_K_LOWTHRESH = 1,
        ADC_V_LOWTHRESH = 13, // determines location in packed byte
        ADC_M_LOWTHRESH = ((1 << ADC_K_LOWTHRESH) - 1),

        ADC_K_HIGHTHRESH = 1,
        ADC_V_HIGHTHRESH = 14,  // determines location in packed byte
        ADC_M_HIGHTHRESH = ((1 << ADC_K_HIGHTHRESH) - 1)
    };



    /// Acd ID
    idents::AcdId       m_ID;
    /// pulse height
    unsigned short       m_pulseHeight;
    /// nominal Acd veto signal
    bool                 m_veto;
    /// 1 bit Low threshold discriminator - enables the PHA
    bool                 m_lowThreshold;
    /// 1 bit High threshold discriminator - used for calibration of the CAL
    bool                 m_highThreshold;
    /// packed PHA word, containing status bits and PHA for PDS
    unsigned short       m_packedPHA;
};


//! Definition of all container types of AcdDigi
template <class TYPE> class ObjectVector;
typedef ObjectVector<AcdDigi>     AcdDigiVector;

template <class TYPE> class ObjectList;
typedef ObjectList<AcdDigi> AcdDigiList;

/// Serialize the object for writing
inline StreamBuffer& AcdDigi::serialize( StreamBuffer& s ) const
{
  ContainedObject::serialize(s);
  return s
    << m_ID
    << m_packedPHA;
}


/// Serialize the object for reading
inline StreamBuffer& AcdDigi::serialize( StreamBuffer& s )
{
  ContainedObject::serialize(s);
  
  unsigned int id;

  s >> id
    >> m_packedPHA;

   idents::AcdId persId(id);
   m_ID = persId;

  // unpack the bits from m_packedPHA
  m_veto =  (m_packedPHA >> ADC_V_VETO) & ADC_M_VETO;
  m_lowThreshold =  (m_packedPHA >> ADC_V_LOWTHRESH) & ADC_M_LOWTHRESH;
  m_highThreshold =  (m_packedPHA >> ADC_V_HIGHTHRESH) & ADC_M_HIGHTHRESH;

  return s;
}


/// Fill the ASCII output stream
inline std::ostream& AcdDigi::fillStream( std::ostream& s ) const
{
  return s
    << "    base class AcdDigi :"
    << "\n        ID = ( "
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_ID << ", "
    << "\n        pulse height      = "
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_pulseHeight << ", "
    << "\n        veto              = "
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_veto   << " )"
    << "\n        Low Threshold     = "
    << m_lowThreshold << " )"
    << "\n        High Threshold     = "
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_highThreshold << " )"
    << "\n        packed PHA        = "
    << EventFloatFormat( EventFormat::width, EventFormat::precision )
    << m_packedPHA << " )";
}




#endif    // GlastEvent_AcdDigi_H

