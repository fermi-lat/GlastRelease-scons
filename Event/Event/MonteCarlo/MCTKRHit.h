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


//------------------------------------------------------------------------------
//
// ClassName:   MCTKRHit
//  
// Description: Essential information of the MCTKRHit
//
//              It contains:
//                  - id
//                  - energy
//                  - noise
//
//
//------------------------------------------------------------------------------

/*!
Essential information of the MCTKRHit.

  It contains:
  - id
  - energy
  - noise
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
    
    //! energy deposited
    float energy () const              { return m_energy; }
    void setEnergy (float value)        { m_energy = value; }

    /// noise
    int noise () const           { return m_noise; }
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
