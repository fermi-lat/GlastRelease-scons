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


//------------------------------------------------------------------------------
//
// ClassName:   MCACDHit
//  
// Description: Essential information of the MCACDHit
//
//              It contains:
//                  - id
//                  - energy
//
//
//------------------------------------------------------------------------------

/*!
Essential information of the MCACDHit.

  It contains:
  - id
  - energy
*/

// declare the GAUDI class-id for this class
//const static CLID   CLID_MCACDHit = 205;


class MCACDHit : virtual public ContainedObject  {  
    
public:
    /// Constructors
    MCACDHit() { }
    
    /// Destructor
    virtual ~MCACDHit() { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return MCACDHit::classID(); }
    static const CLID& classID()       { return CLID_MCACDHit; }
    
    /// Tile id number
    unsigned int id () const           { return m_id; }
    void setId (long value)            { m_id = value; }
    
    /// energy deposited (GeV)
    float energy () const              { return m_energy; }
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


// Definition of all container types of MCACDHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<MCACDHit>     MCACDHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<MCACDHit>       MCACDHitList;


#endif    // MCACDHit_H
