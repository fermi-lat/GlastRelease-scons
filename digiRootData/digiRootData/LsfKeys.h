/** 
 *
 * @brief Classes encapsulating translated file keys from event data 
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef ROOT_LsfKeys_H
#define ROOT_LsfKeys_H

#include "Riostream.h"
#include "TObject.h"

#include <vector>
//#include <iostream>
//#include <iomanip>

#include "enums/Lsf.h"


// forward declare the subclasses
class LpaKeys;
class LciKeys;

// abstract base keys class
class LsfKeys : public TObject {

  public:
    LsfKeys() : m_LATC_master(0), m_LATC_ignore(0) {};

    LsfKeys( UInt_t master, UInt_t ignore ) 
      : m_LATC_master( master ), m_LATC_ignore( ignore ) {};

    LsfKeys( const LsfKeys& other )
      : TObject(other), m_LATC_master( other.LATC_master() ), m_LATC_ignore( other.LATC_ignore() ) {}

    virtual ~LsfKeys() {};

    void Print(Option_t* /*option=""*/) const;
    void Fake(Int_t ievent, Float_t randNum);
    Bool_t CompareInRange( const LsfKeys&, const std::string & name = "" ) const ; // for tests

    const LsfKeys& operator=( const LsfKeys& other ) {
      if ( &other != this ) {
	m_LATC_master = other.LATC_master();
	m_LATC_ignore = other.LATC_ignore();
      }
      return *this;
    };

    UInt_t LATC_master() const { return m_LATC_master; };
    UInt_t LATC_ignore() const { return m_LATC_ignore; };

    void setLATC_master( UInt_t value ) { m_LATC_master = value; };
    void setLATC_ignore( UInt_t value ) { m_LATC_ignore = value; };


    virtual LsfKeys* clone() const = 0;
    virtual enums::Lsf::KeysType type() const = 0;
    virtual const LpaKeys* castToLpaKeys() const { return 0; };
    virtual const LciKeys* castToLciKeys() const { return 0; };

  private:
    UInt_t m_LATC_master;
    UInt_t m_LATC_ignore;

    friend class LpaKeys;
    friend class LciKeys;

    ClassDef(LsfKeys,1) // Base class for LAT keys
  };

  // translated keys from LPA data
  class LpaKeys : public LsfKeys {

  public:
    LpaKeys() : LsfKeys(), m_sbs(0), m_lpa_db(0) {};

    LpaKeys( UInt_t master, UInt_t ignore, UInt_t sbs, UInt_t lpa_db)
      : LsfKeys( master, ignore ), m_sbs(sbs), m_lpa_db(lpa_db) {};

    LpaKeys( const LpaKeys& other ) : LsfKeys( other ), m_sbs( other.sbs() ),
                                      m_lpa_db(other.lpa_db()) {};

    virtual ~LpaKeys() { };

    void Print(Option_t* option="" ) const;
    void Fake(Int_t ievent, Float_t randNum);
    Bool_t CompareInRange( const LpaKeys&, const std::string & name = "" ) const ; // for tests

    const LpaKeys& operator=( const LpaKeys& other ) {
      if ( &other != this ) {
	( dynamic_cast< LsfKeys& >( *this ) ) = other;
      }
      return *this;
    }

    void setSbs(UInt_t val) { m_sbs = val; }
    UInt_t sbs() const { return m_sbs; }

    void setLpa_db(UInt_t val) { m_lpa_db = val ; }
    UInt_t lpa_db() const { return m_lpa_db; }

    virtual LsfKeys*  clone() const { return new LpaKeys( *this ); };
    virtual const LpaKeys* castToLpaKeys() const { return this; };
    virtual enums::Lsf::KeysType type() const { return enums::Lsf::LpaKeys; };

  private:
    std::vector<UInt_t> m_CDM_keys;
    UInt_t m_sbs; /// FMX key of the secondary boot script
    UInt_t m_lpa_db; /// same as LPA_Info::lpaDbKey
    ClassDef(LpaKeys,2) // LPA Key
  };

  // keys from LCI data
  class LciKeys : public LsfKeys {

  public:
    LciKeys() : LsfKeys(), m_LCI_script( 0 ) {};

    LciKeys( UInt_t master, UInt_t ignore, UInt_t script )
      : LsfKeys( master, ignore ), m_LCI_script( script ) {};

    LciKeys( const LciKeys& other ) : LsfKeys( other ), m_LCI_script( other.LCI_script() ) {};

    virtual ~LciKeys() {};

    void Print(Option_t* option="" ) const;
    void Fake(Int_t ievent, Float_t randNum);
    Bool_t CompareInRange( const LciKeys&, const std::string & name = "" ) const ; // for tests

    const LciKeys& operator=( const LciKeys& other ) {
      if ( &other != this ) {
	( dynamic_cast< LsfKeys& >( *this ) ) = other;
	m_LCI_script = other.LCI_script();
      }
      return *this;
    };

    UInt_t LCI_script() const { return m_LCI_script; };

    void setLCI_script( UInt_t value ) { m_LCI_script = value; };


    virtual LsfKeys*  clone() const { return new LciKeys( *this ); };
    virtual enums::Lsf::KeysType type() const { return enums::Lsf::LciKeys; };
    virtual const LciKeys* castToLciKeys() const { return this; };

  private:
    UInt_t m_LCI_script;
    ClassDef(LciKeys,1) // LCI Keys
  };


#endif
