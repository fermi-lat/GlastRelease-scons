
#ifndef IrfAcdHit_H
#define IrfAcdHit_H 1

// Include files
#include <iostream>
#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/StreamBuffer.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "GlastEvent/TopLevel/Definitions.h"

// Externals 
extern const CLID& CLID_IrfAcdHit;

/*! \class IrfAcdHit
\brief IRF data concerning the ACD 

  It contains:
  - ACD tile id
  - energy deposited (GeV)
*/


class IrfAcdHit : virtual public ContainedObject  {  
    
public:
    /// Constructor
    IrfAcdHit() { }
    
    /// Destructor
    virtual ~IrfAcdHit() { }
    
    //! Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return IrfAcdHit::classID(); }
    static const CLID& classID()       { return CLID_IrfAcdHit; }
    
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

inline StreamBuffer& IrfAcdHit::serialize( StreamBuffer& s ) const                 {
    ContainedObject::serialize(s);
    return s
	<< m_id
	<< m_energy
	;
}


inline StreamBuffer& IrfAcdHit::serialize( StreamBuffer& s )                       {
    ContainedObject::serialize(s);
    return s
	>> m_id
	>> m_energy
	;
}


inline std::ostream& IrfAcdHit::fillStream( std::ostream& s ) const                {
    return s
	<< "class IrfAcdHit :"
	<< "\n    Energy    = "
	//  << LHCbEventFloatFormat( LHCbEvent::width, LHCbEvent::precision )
	<< m_energy
	<< "\n    Cell ID   = " << m_id;
}


//! Definition of all container types of IrfAcdHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<IrfAcdHit>     IrfAcdHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<IrfAcdHit>       IrfAcdHitList;


#endif    // IrfAcdHit_H
