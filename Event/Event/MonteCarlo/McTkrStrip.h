#ifndef McTkrStrip_H
#define McTkrStrip_H 1

//include files

#include <iostream>
#include <vector>

//#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"

#include "Event/TopLevel/Definitions.h"
#include "GaudiKernel/ObjectVector.h"
#include "idents/VolumeIdentifier.h"

/*!
* \class McTkrStrip
* \author T. Burnett
* \brief Represent energy deposited (or charge) per Si Strip
*
*  Intermediate class for TkrDigi process.
*/

extern const CLID& CLID_McTkrStrip;

namespace Event {
  class McTkrStrip : virtual public ContainedObject {
    
  public:
    //! typedefs
    
    //! Constructors
    //! Null constructor
    //McTkrStrip() {};
    
    //! constructor with plane id.
    McTkrStrip(idents::VolumeIdentifier id, unsigned int strip,  double e=0 )
        : m_planeId(id)
        , m_strip(strip)
        , m_energy(e) {};
    //! Destructor
    virtual ~McTkrStrip() {};
        
    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID() const   { return McTkrStrip::classID(); }
    static const CLID& classID()       { return CLID_McTkrStrip; }
    
    //! access to the energy
    double getEnergy()const { return m_energy; }
    double energy()const { return m_energy; }
    
    //! increase the energy
    void operator+=(double de){m_energy += de;}

    //! access to the id
    idents::VolumeIdentifier getId()const { return m_planeId; }

    //! access to the strip number
    unsigned int getStripNumber()const { return m_strip; }

    //! Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    //! Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    //! Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;        
        
  private:
    //! the plane that this belongs to
    idents::VolumeIdentifier m_planeId;
    
    //! strip number
    unsigned int m_strip;
    
    //! total energy deposited
    double m_energy;
  };
    
    
  //! Serialize the object for writing
  inline StreamBuffer& McTkrStrip::serialize( StreamBuffer& s ) const {
    ContainedObject::serialize(s);  
#if 0
    s   << m_planeId
        << m_strip
        << m_energy;
#endif 
    return s;
  }
  
  //! Serialize the object for reading
  inline StreamBuffer& McTkrStrip::serialize( StreamBuffer& s )       {
    ContainedObject::serialize(s);
#if 0
    s   >> m_planeId
        >> m_strip
        >> tower
        >> m_energy;
#endif
    return s;
  }
  
  //! Fill the ASCII output stream
  
  inline std::ostream& McTkrStrip::fillStream( std::ostream& s ) const {
    s   << "class TkrStrip :" << std::endl
        << "Plane Id: " << m_planeId.name() 
        << "   strip: " << m_strip
        << "  energy: " << m_energy 
        << std::endl;
    
    return s;
  }
  
  //! Definition of all container type of TkrStrip
  
  typedef ObjectVector<McTkrStrip> McTkrStripCol;
} //Namespace Event
#endif
