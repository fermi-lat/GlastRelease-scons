// $Header$
// Author: H. Arrighi

#ifndef ACDhit_H
#define ACDhit_H 1

// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "GlastEvent/TopLevel/Definitions.h"

// Externals 
extern const CLID& CLID_ACDhit;


//------------------------------------------------------------------------------
//
// ClassName:   ACDhit
//  
// Description: Essential information of the ACDhit
//
//              It contains:
//                  - id
//                  - energy
//
//
//------------------------------------------------------------------------------

/*!
Essential information of the ACDhit.

  It contains:
  - id
  - energy
*/

// declare the GAUDI class-id for this class
//const static CLID   CLID_ACDhit = 205;


class ACDhit : virtual public ContainedObject  {  
    
public:
    /// Constructors
    ACDhit() { }
    
    /// Destructor
    virtual ~ACDhit() { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return ACDhit::classID(); }
    static const CLID& classID()       { return CLID_ACDhit; }
    
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

inline StreamBuffer& ACDhit::serialize( StreamBuffer& s ) const                 {
    ContainedObject::serialize(s);
    return s
	<< m_id
	<< m_energy
	;
}


inline StreamBuffer& ACDhit::serialize( StreamBuffer& s )                       {
    ContainedObject::serialize(s);
    return s
	>> m_id
	>> m_energy
	;
}


inline std::ostream& ACDhit::fillStream( std::ostream& s ) const                {
    return s
	<< "class ACDhit :"
	<< "\n    Energy    = "
	//  << LHCbEventFloatFormat( LHCbEvent::width, LHCbEvent::precision )
	<< m_energy
	<< "\n    Cell ID   = " << m_id;
}


// Definition of all container types of ACDhit
template <class TYPE> class ObjectVector;
typedef ObjectVector<ACDhit>     ACDhitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<ACDhit>       ACDhitList;


#endif    // ACDhit_H
