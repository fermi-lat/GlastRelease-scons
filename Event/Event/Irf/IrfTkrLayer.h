#ifndef IrfTkrLayer_H
#define IrfTkrLayer_H 1

// Include files
#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GlastEvent/TopLevel/Definitions.h"

#include "GlastEvent/Irf/IrfTkrHit.h"

// Externals 
extern const CLID& CLID_IrfTkrLayer;


//------------------------------------------------------------------------------
//
// ClassName:   IrfTkrLayer
//  
// Description: Essential information of the IrfTkrLayer
//
//              It contains:
//                  - id
//                  - max energy energy
//                  - vector of hit strips
//
//
//------------------------------------------------------------------------------

/*!
Essential information of the IrfTkrLayer.

  It contains:
  - id
  - max energy
  - vector of hit strips
*/


class IrfTkrLayer : virtual public ContainedObject  {  
    
public:
    /// Constructors
    IrfTkrLayer() { }
    
    /// Destructor
    virtual ~IrfTkrLayer() { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return IrfTkrLayer::classID(); }
    static const CLID& classID()       { return CLID_IrfTkrLayer; }
    
    //! SiLayer id number
    unsigned int id () const           { return m_id; }
    void setId (long value)            { m_id = value; }
    
    //! max energy deposited in all hit strips for this SiLayer
    float MaxEnergy () const              { return m_energy; }
    void setMaxEnergy (float value)        { m_energy = value; }

    //! clear hitVector
    void clearHits () {
        for (IrfTkrHitVector::iterator it = m_hitVector.begin(); it != m_hitVector.end(); it++) {
            m_hitVector.erase(it);
        }
    }

    //! add a IrfTkrHit to the list of hit strips for this SiLayer
   void addHit (IrfTkrHit *hit) { m_hitVector.push_back(hit); }

    //! retrieve the vector of hit strips
    IrfTkrHitVector* getHits() { return &m_hitVector; }

    //! used to convert container of objects
    virtual StreamBuffer& serialize( StreamBuffer& s );
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
private:
    unsigned int         m_id;        // ID of layer
    float                m_energy;      // max energy in layer
    IrfTkrHitVector      m_hitVector;  // list of hit strips
};

inline StreamBuffer& IrfTkrLayer::serialize( StreamBuffer& s ) const                 {
    s = ContainedObject::serialize(s);
    s
        << m_id
	<< m_energy
        << m_hitVector.size();
    for (IrfTkrHitVector::const_iterator hit = m_hitVector.begin(); hit != m_hitVector.end(); hit++) {
         s = (*hit)->serialize(s);
    }
    return s;
}


inline StreamBuffer& IrfTkrLayer::serialize( StreamBuffer& s )                       {
    s = ContainedObject::serialize(s);
    int nSize;
    s
	>> m_id
	>> m_energy
        >> nSize;
    for (int iHit = 0; iHit < nSize; iHit++){
        IrfTkrHit* curHit = new IrfTkrHit();
        s = curHit->serialize(s);
        this->addHit(curHit);
    }

    return s;
}


inline std::ostream& IrfTkrLayer::fillStream( std::ostream& s ) const                {
    return s
	<< "class IrfTkrLayer :"
	<< "\n    Max Energy    = "
	<< m_energy
	<< "\n    Si Layer id   = " << m_id
        << "\n    hit strips    = " << m_hitVector;
}

// Definition of all container types of IrfTkrLayer
template <class TYPE> class ObjectVector;
typedef ObjectVector<IrfTkrLayer>     IrfTkrLayerVector;
template <class TYPE> class ObjectList;
typedef ObjectList<IrfTkrLayer>       IrfTkrLayerList;


#endif    // IrfTkrLayer_H
