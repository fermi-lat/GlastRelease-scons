#ifndef GlastEvent_AcdDigi_H
#define GlastEvent_AcdDigi_H 1


// Include files
#include <iostream>
#include "idents/ScintillatorId.h"
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/ContainedObject.h"

#include "GlastEvent/TopLevel/Definitions.h"

/*!
//------------------------------------------------------------------------------
//
// \class   AcdDigi        
//  
// \brief Digitizations for ACD                                
//              
// Author:  Richard Dubois 
//
// Note the usage of ScintillatorId, this will be replaced by
// the new AcdId class under development.
//
// There are no set methods in this class, users are expected to fill
// the data members through the constructor.
//------------------------------------------------------------------------------
 */

extern const CLID& CLID_AcdDigi;


class AcdDigi : virtual public ContainedObject  { 

public:
    AcdDigi() {};

    AcdDigi(idents::ScintillatorId id, int pulseHeight, bool veto, bool lowThresh, 
        bool highThresh)       
        : m_ID(id),
        m_pulseHeight(pulseHeight),
        m_veto(veto),
        m_lowThreshold(lowThresh),
        m_highThreshold(highThresh)
    {};

    /// Destructor
    virtual ~AcdDigi() { };

    //! Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return AcdDigi::classID(); }
    static const CLID& classID()       { return CLID_AcdDigi; }

    /// Retrieve digi identifier
    const idents::ScintillatorId ID() const;

    /// Retrieve pulse height
    const int PulseHeight() const;

    const bool Veto() const;

    const bool LowThreshold() const;

    const bool HighThreshold() const;


    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;



private:


    enum {
        ADC_K_VETO = 1,
        ADC_V_VETO = 12,
        ADC_M_VETO = ((1 << ADC_K_VETO) - 1),

        ADC_K_LOWTHRESH = 1,
        ADC_V_LOWTHRESH = 13,
        ADC_M_LOWTHRESH = ((1 << ADC_K_LOWTHRESH) - 1),

        ADC_K_HIGHTHRESH = 1,
        ADC_V_HIGHTHRESH = 14,
        ADC_M_HIGHTHRESH = ((1 << ADC_K_HIGHTHRESH) - 1)
    };



    /// Acd ID
    idents::ScintillatorId       m_ID;
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
 // s >> m_ID  need to add for ScintillatorId
   s >> m_packedPHA;

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
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_ID << ", "
    << "\n        pulse height      = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_pulseHeight << ", "
    << "\n        veto              = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_veto   << " )"
    << "\n        Low Threshold     = "
    << m_lowThreshold << " )"
    << "\n        High Threshold     = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_highThreshold << " )"
    << "\n        packed PHA        = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_packedPHA << " )";
}




#endif    // GlastEvent_AcdDigi_H

