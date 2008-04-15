// -*- Mode: C++; -*-
/** 
 * @namespace eventFile
 *
 * @brief Classes encapsulating translated file keys from event data 
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef EVENTFILE_LSE_KEYS_HH
#define EVENTFILE_LSE_KEYS_HH

#include <vector>

namespace eventFile {

  struct LSE_Keys {
    typedef enum _KeysType {
      NONE = -1,
      LPA  = 0,
      LCI,
      NumLSEKeysTypes,
    } KeysType;

    LSE_Keys( unsigned master, unsigned ignore)
      : LATC_master( master ), LATC_ignore( ignore ) {};
    virtual ~LSE_Keys() {};
    virtual void dump( const char* pre, const char* post ) const;
    unsigned LATC_master;
    unsigned LATC_ignore;
  };

  struct LPA_Keys : public LSE_Keys {
    LPA_Keys( unsigned master=0xFFFFFFFF, unsigned ignore=0xFFFFFFFF, 
	      unsigned sbs=0xFFFFFFFF, unsigned lpadb=0xFFFFFFFF )
      : LSE_Keys( master, ignore ), SBS( sbs ), LPA_db( lpadb ) {};
    void dump( const char* pre, const char* post ) const;
    unsigned SBS;
    unsigned LPA_db;
  };

  struct LCI_Keys : public LSE_Keys {
    LCI_Keys( unsigned master=0xFFFFFFFF, unsigned ignore=0xFFFFFFFF, unsigned script=0xFFFFFFFF )
      : LSE_Keys( master, ignore ), LCI_script( script ) {};
    void dump( const char* pre, const char* post ) const;
    unsigned LCI_script;
  };
  
};

#endif // EVENTFILE_LSE_KEYS_HH

