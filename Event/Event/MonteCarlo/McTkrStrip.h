#ifndef McTkrStrip_H
#define McTkrStrip_H 1

//include files

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"

#include "Event/TopLevel/Definitions.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"

#include "idents/TowerId.h"
#include "idents/GlastAxis.h"

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
    McTkrStrip(idents::VolumeIdentifier id(m_planeId), 
               unsigned int strip(m_strip)) {};
    //! Destructor
    virtual ~McTkrStrip() {};
        
    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID() const   { return TkrStrip::classID(); }
    static const CLID& classID()       { return CLID_TkrStrip; }
    
    //! access to the energy
    double getEnergy()const { return m_energy; }
    
    //! increase the energy
    void operator+=(double de){m_energy += de);

    //! access to the id
    double getId()const { return m_planeId; }

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
    unsigned in m_strip;
    
    //! total energy deposited
    double m_energy;
  };
    
    
  //! Serialize the object for writing
  inline StreamBuffer& McTkrStrip::serialize( StreamBuffer& s ) const {
    ContainedObject::serialize(s);  
    s   << m_planeId
        << m_strip
        << m_energy;
    
    return s;
  }
  
  //! Serialize the object for reading
  inline StreamBuffer& McTkrStrip::serialize( StreamBuffer& s )       {
    ContainedObject::serialize(s);
    s   >> m_planeId
        >> m_strip
        >> tower
        >> m_energy;
    
    return s;
  }
  
  //! Fill the ASCII output stream
  
  inline std::ostream& McTkrStrip::fillStream( std::ostream& s ) const {
    int j;
    int size = m_hits.size();
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
