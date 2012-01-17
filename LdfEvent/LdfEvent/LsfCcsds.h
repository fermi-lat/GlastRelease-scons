#ifndef LDFEVENT_CCSDS_H
#define LDFEVENT_CCSDS_H 1

#include <cstdio>
#include <iostream>
#include "GaudiKernel/StreamBuffer.h"

#include "Event/TopLevel/Definitions.h"

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "lsfData/LsfCcsds.h"

/** @class LsfCcsds 
*
* $Header$
*/

static const CLID& CLID_LsfCcsds = InterfaceID("LsfCcsds", 1, 0);

namespace LsfEvent {
  
  class LsfCcsds : public DataObject {
    
  public:
    
    /// Default c'tor.  Assigns sentinel or null values to all fields
    LsfCcsds():DataObject() { clear(); }

    /// Standard c'tor.  Takes input values for all fields
    LsfCcsds( int scid, int apid, double utc)
      :DataObject(),m_scid(scid),m_apid(apid),m_utc(utc) {
    }

    /// Copy c'tor.  Just copy all values.  
    LsfCcsds( const LsfCcsds& other )
      :DataObject(other),
       m_scid(other.m_scid),
       m_apid(other.m_apid),
       m_utc(other.m_utc) {
    }

    /// Copy c'tor.  Just copy all values.  
    LsfCcsds( const lsfData::LsfCcsds& other )
      :DataObject(),
       m_scid(other.getScid()),
       m_apid(other.getApid()),
       m_utc(other.getUtc()) {
     }
    
    /// D'tor
    virtual ~LsfCcsds(){
    }
    
    void clear() {
        m_scid = 0;
        m_apid = 0;
        m_utc = 0.0;
    }

    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return LsfCcsds::classID(); }
    static const CLID& classID()       { return CLID_LsfCcsds; }
    
    int scid() const { return m_scid; };
    int apid() const { return m_apid; };
    double utc() const { return m_utc; };

    /// set everything at once
    inline void initialize(int scid, int apid, double utc) {
      m_scid = scid;
      m_apid = apid;
      m_utc = utc;
    }

    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const LsfCcsds& obj ) {
       return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const {
      s << "class LsfCcsds : \n"
	<< EventField( EventFormat::field12 )
	<< m_scid << "\n"
	<< m_apid << "\n"
	<< m_utc << "\n";
      return s;
    }
    
  private:
    
      int m_scid;
      int m_apid;
      double m_utc;

  };

}


#endif    // LSF_CCSDS
