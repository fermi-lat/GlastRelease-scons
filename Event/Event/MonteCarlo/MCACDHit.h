// $Header$
// Author: H. Arrighi

#ifndef MCACDHit_H
#define MCACDHit_H 1

// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "GlastEvent/TopLevel/Definitions.h"

// Externals 
extern const CLID& CLID_MCACDHit;

/*! \class MCACDHit
\brief MC data concerning the ACD - currently filled with IRF data

  It contains:
  - ACD tile id
  - energy deposited (GeV)
*/


class MCACDHit : virtual public ContainedObject  {  
    
public:
    /// Constructor
    MCACDHit() { }
    
    /// Destructor
    virtual ~MCACDHit() { }
    
    //! Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return MCACDHit::classID(); }
    static const CLID& classID()       { return CLID_MCACDHit; }
    
    //! Tile id number
    unsigned int id () const           { return m_id; }
    void setId (long value)            { m_id = value; }
    
    //! get energy (GeV) deposited in this ACD tile
    float energy () const              { return m_energy; }
    //! set the energy (GeV) deposited in this ACD tile
    void setEnergy (float value)        { m_energy = value; }
    
    //! used to convert container of objects
    virtual StreamBuffer& serialize( StreamBuffer& s );
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
private:
    unsigned int         m_id;
    float                m_energy;
};

inline StreamBuffer& MCACDHit::serialize( StreamBuffer& s ) const                 {
    ContainedObject::serialize(s);
    return s
	<< m_id
	<< m_energy
	;
}


inline StreamBuffer& MCACDHit::serialize( StreamBuffer& s )                       {
    ContainedObject::serialize(s);
    return s
	>> m_id
	>> m_energy
	;
}


inline std::ostream& MCACDHit::fillStream( std::ostream& s ) const                {
    return s
	<< "class MCACDHit :"
	<< "\n    Energy    = "
	//  << LHCbEventFloatFormat( LHCbEvent::width, LHCbEvent::precision )
	<< m_energy
	<< "\n    Cell ID   = " << m_id;
}


//! Definition of all container types of MCACDHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<MCACDHit>     MCACDHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<MCACDHit>       MCACDHitList;


#endif    // MCACDHit_H
