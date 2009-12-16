#ifndef Event_AcdReconV2_H
#define Event_AcdReconV2_H 1

#include <iostream>
#include "idents/AcdId.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/DataObject.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/Recon/AcdRecon/AcdTkrAssoc.h"
#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Recon/AcdRecon/AcdEventTopology.h"

#include <vector>

#include "GaudiKernel/IInterface.h"

static const CLID& CLID_AcdReconV2 = InterfaceID("AcdReconV2", 1, 0);

/** @class Event::AcdReconV2        
* @brief Reconstruction data for ACD
*
* The reconstruction data consists of:
*
*                                 
* @author Heather Kelly
* $Header$          
*/

namespace Event {
    
  class AcdReconV2 : virtual public DataObject  { 
    
  public:
    
    static AcdReconV2* s_theAcdReconV2Ptr;
    
  public:

    AcdReconV2(){};
    
    virtual ~AcdReconV2(){};
    
    void init();
    
    void clear();

    //! Retrieve reference to class definition structure
    virtual const CLID& clID() const    { return AcdReconV2::classID(); }
    static const  CLID& classID()       { return CLID_AcdReconV2; }
    

    // Access to the data
    inline AcdEventTopology& getEventTopology() { return m_topology; }
    inline const AcdEventTopology& getEventTopology() const { return m_topology; }

    inline AcdHitCol& getHitCol() { return m_acdHits; }
    inline const AcdHitCol& getHitCol() const { return m_acdHits; }
    
    inline AcdTkrAssocCol& getTkrAssocCol() { return m_acdTkrAssocs; }
    inline const AcdTkrAssocCol& getTkrAssocCol() const { return m_acdTkrAssocs; }
    
    
    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    
    
    friend std::ostream& operator << (std::ostream& s, const AcdReconV2& obj) {
      return obj.fillStream(s);
    };
    
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
  private:

    /// Event Topology
    AcdEventTopology m_topology;
        
    /// the vector of calibrated ACD hits
    AcdHitCol m_acdHits;	
    
    /// Track associations
    AcdTkrAssocCol m_acdTkrAssocs;

  };

  inline void AcdReconV2::init() {
    clear();
  }

  inline void AcdReconV2::clear() {
    m_topology.ini();
    m_acdHits.clear();
    m_acdTkrAssocs.clear();    
  }

  
  
  /// Serialize the object for writing
  inline StreamBuffer& AcdReconV2::serialize( StreamBuffer& s ) const {
    DataObject::serialize(s);
    return s;
  }
    
    
    /// Serialize the object for reading
  inline StreamBuffer& AcdReconV2::serialize( StreamBuffer& s ) {
    DataObject::serialize(s);
    return s;
  }
    
    
  /// Fill the ASCII output stream
  inline std::ostream& AcdReconV2::fillStream( std::ostream& s ) const {
    return s;
  }
        
} // namespace Event

#endif    // Event_AcdReconV2_H
