#ifndef MCTKRHit_H
#define MCTKRHit_H 1

// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "GlastEvent/TopLevel/Definitions.h"

// Externals 
extern const CLID& CLID_MCTKRHit;

/*! \class MCTKRHit
\brief MC class for the TKR  - currently filled with IRF data

  It contains:
  - strip id
  - energy deposited
  - noise, denotes whether (1) or not (0) this hit was caused by noise
*/


class MCTKRHit : virtual public ContainedObject  {  
    
public:
    /// Constructors
    MCTKRHit() { 
        m_id = 0;
        m_energy = 0.f;
        m_noise = 0;
    }
    
    /// Destructor
    virtual ~MCTKRHit() { }
    
    //! Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return MCTKRHit::classID(); }
    static const CLID& classID()       { return CLID_MCTKRHit; }
    
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

inline StreamBuffer& MCTKRHit::serialize( StreamBuffer& s ) const                 {
    ContainedObject::serialize(s);
    return s
	<< m_id
	<< m_energy
        << m_noise;
}


inline StreamBuffer& MCTKRHit::serialize( StreamBuffer& s )                       {
    ContainedObject::serialize(s);
    return s
	>> m_id
	>> m_energy
        >> m_noise;
}


inline std::ostream& MCTKRHit::fillStream( std::ostream& s ) const                {
    return s
	<< "class MCTKRHit :"
	<< "\n    Energy    = "
	<< m_energy
	<< "\n    Strip id   = " << m_id
        << "\n    Noise = " << m_noise;
}


// Definition of all container types of MCTKRHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<MCTKRHit>     MCTKRHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<MCTKRHit>       MCTKRHitList;


#endif    // MCTKRHit_H
