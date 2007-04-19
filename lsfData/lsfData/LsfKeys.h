// -*- Mode: C++; -*-
/** 
 * @namespace lsfData
 *
 * @brief Classes encapsulating translated file keys from event data 
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef lsfData_LsfKeys_H
#define lsfData_LsfKeys_H

#include <vector>
#include <iostream>
#include <iomanip>

#include "enums/Lsf.h"

namespace lsfData {

  // forward declare the subclasses
  class LpaKeys;
  class LciKeys;

  // abstract base keys class
  class LsfKeys {
    unsigned int m_LATC_master;
    unsigned int m_LATC_ignore;

  public:
    LsfKeys() : m_LATC_master(0), m_LATC_ignore(0) {};

    LsfKeys( unsigned int master, unsigned int ignore ) 
      : m_LATC_master( master ), m_LATC_ignore( ignore ) {};

    LsfKeys( const LsfKeys& other )
      : m_LATC_master( other.LATC_master() ), m_LATC_ignore( other.LATC_ignore() ) {}

    virtual ~LsfKeys() {};

    const LsfKeys& operator=( const LsfKeys& other ) {
      if ( &other != this ) {
	m_LATC_master = other.LATC_master();
	m_LATC_ignore = other.LATC_ignore();
      }
      return *this;
    };

    unsigned int LATC_master() const { return m_LATC_master; };
    unsigned int LATC_ignore() const { return m_LATC_ignore; };

    void setLATC_master( unsigned int value ) { m_LATC_master = value; };
    void setLATC_ignore( unsigned int value ) { m_LATC_ignore = value; };

    virtual void print( const std::string& str="" ) const {
      std::cout << str << " LATC_master = 0x" << std::setw(8) << std::setfill('0') << std::hex << m_LATC_master << std::endl;
      std::cout << str << " LATC_ignore = 0x" << std::setw(8) << std::setfill('0') << std::hex << m_LATC_ignore << std::endl;
    }

    virtual LsfKeys* clone() const = 0;
    virtual enums::Lsf::KeysType type() const = 0;
    virtual const LpaKeys* castToLpaKeys() const { return 0; };
    virtual const LciKeys* castToLciKeys() const { return 0; };

  };

  // translated keys from LPA data
  class LpaKeys : public LsfKeys {
    std::vector< unsigned int > m_LPA_DB;

  public:
    LpaKeys() : LsfKeys(), m_LPA_DB( std::vector<unsigned int>() ) {};

    LpaKeys( unsigned int master, unsigned int ignore, const std::vector< unsigned int >& dbs )
      : LsfKeys( master, ignore ), m_LPA_DB( dbs ) {};

    LpaKeys( const LpaKeys& other ) : LsfKeys( other ), m_LPA_DB( other.LPA_DB() ) {};

    virtual ~LpaKeys() {};

    const LpaKeys& operator=( const LpaKeys& other ) {
      if ( &other != this ) {
	( dynamic_cast< LsfKeys& >( *this ) ) = other;
	m_LPA_DB.clear();
	m_LPA_DB.insert( m_LPA_DB.begin(), other.LPA_DB().begin(), other.LPA_DB().end() );
      }
      return *this;
    }

    const std::vector< unsigned int >& LPA_DB() const { return m_LPA_DB; };

    void setLPA_DB( const std::vector< unsigned int >& value ) {
      m_LPA_DB.clear();
      m_LPA_DB.insert( m_LPA_DB.begin(), value.begin(), value.end() );
    };

    virtual void print( const std::string& str="" ) const {
      LsfKeys::print( str );
      std::vector<unsigned int>::const_iterator itr = m_LPA_DB.begin();
      for ( int i=0; itr != m_LPA_DB.end(); ++itr, ++i ) {
	std::cout << str << " LPA_DB[" << i << "] = 0x" << std::setw(8) << std::setfill('0') << std::hex << *itr << std::endl;
      }
    }
    
    virtual LsfKeys*                   clone() const { return new LpaKeys( *this ); };
    virtual const LpaKeys*  castToLpaKeys() const { return this; };
    virtual enums::Lsf::KeysType        type() const { return enums::Lsf::LpaKeys; };

  };

  // keys from LCI data
  class LciKeys : public LsfKeys {
    unsigned int m_LCI_script;

  public:
    LciKeys() : LsfKeys(), m_LCI_script( 0 ) {};

    LciKeys( unsigned int master, unsigned int ignore, unsigned int script )
      : LsfKeys( master, ignore ), m_LCI_script( script ) {};

    LciKeys( const LciKeys& other ) : LsfKeys( other ), m_LCI_script( other.LCI_script() ) {};

    virtual ~LciKeys() {};

    const LciKeys& operator=( const LciKeys& other ) {
      if ( &other != this ) {
	( dynamic_cast< LsfKeys& >( *this ) ) = other;
	m_LCI_script = other.LCI_script();
      }
      return *this;
    };

    unsigned int LCI_script() const { return m_LCI_script; };

    void setLCI_script( unsigned int value ) { m_LCI_script = value; };

    virtual void print( const std::string& str="" ) const {
      LsfKeys::print( str );
      std::cout << str << " LCI_script = 0x" << std::setw(8) << std::setfill('0') << std::hex << m_LCI_script << std::endl;

    }

    virtual LsfKeys*                   clone() const { return new LciKeys( *this ); };
    virtual enums::Lsf::KeysType        type() const { return enums::Lsf::LciKeys; };
    virtual const LciKeys*  castToLciKeys() const { return this; };

  };

};

#endif
