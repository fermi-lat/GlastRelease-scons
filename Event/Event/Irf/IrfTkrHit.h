#ifndef IrfTkrHit_H
#define IrfTkrHit_H 1

// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "GlastEvent/TopLevel/Definitions.h"

// Externals 
extern const CLID& CLID_IrfTkrHit;

/*! \class IrfTkrHit
\brief Irf class for the TKR

  It contains:
  - strip id
  - energy deposited
  - noise, denotes whether (1) or not (0) this hit was caused by noise
*/


class IrfTkrHit : virtual public ContainedObject  {  
    
public:
    /// Constructors
    IrfTkrHit() { 
        m_id = 0;
        m_energy = 0.f;
        m_noise = 0;
    }
    
    /// Destructor
    virtual ~IrfTkrHit() { }
    
    //! Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return IrfTkrHit::classID(); }
    static const CLID& classID()       { return CLID_IrfTkrHit; }
    
    //! strip id number
    unsigned int id () const           { return m_id; }
    void setId (long value)            { m_id = value; }
    
    //! get energy deposited in this strip
    float energy () const              { return m_energy; }
    //! set the energy deposited in this strip
    void setEnergy (float value)        { m_energy = value; }

    /// get noise value for this strip, where 1 == caused by noise, and 0 == not caused by noise
    int noise () const           { return m_noise; }
    //! set the noise for this strip, where 1 == caused by noise, and 0 == not caused by noise
    void setNoise (int value)   { m_noise = value; }
    
    //! used to convert container of objects
    virtual StreamBuffer& serialize( StreamBuffer& s );
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
private:
    unsigned int         m_id;
    float                m_energy;
    int                  m_noise;
};

inline StreamBuffer& IrfTkrHit::serialize( StreamBuffer& s ) const                 {
    ContainedObject::serialize(s);
    return s
	<< m_id
	<< m_energy
        << m_noise;
}


inline StreamBuffer& IrfTkrHit::serialize( StreamBuffer& s )                       {
    ContainedObject::serialize(s);
    return s
	>> m_id
	>> m_energy
        >> m_noise;
}


inline std::ostream& IrfTkrHit::fillStream( std::ostream& s ) const                {
    return s
	<< "class IrfTkrHit :"
	<< "\n    Energy    = "
	<< m_energy
	<< "\n    Strip id   = " << m_id
        << "\n    Noise = " << m_noise;
}


// Definition of all container types of IrfTkrHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<IrfTkrHit>     IrfTkrHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<IrfTkrHit>       IrfTkrHitList;


#endif    // IrfTkrHit_H
