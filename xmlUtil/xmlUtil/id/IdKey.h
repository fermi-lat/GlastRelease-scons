// $Header$
#ifndef XMLUTIL_IDKEY_H
#define XMLUTIL_IDKEY_H

#include <vector>
#include <string>

#include "xmlUtil/id/NamedId.h"

namespace xmlUtil {

  //*  IdKey is a compressed, sortable version of an identifier which may
  //*  be used as a key for a map.
  class IdKey {
  public:
    //* Construct an IdKey from a "standard" identifier
    IdKey(const Identifier& id);

    /*
    // Construct an IdKey from an identifier whose fields are expressed
    // as strings         NYI
    //    IdKey(const std::vector<std::string>&);  
   */
    ~IdKey();

    //* The comparison function needed when making a map
    class ltkey { 
    public:
      bool operator()(IdKey key1, IdKey key2) const;
    };

  private:
    typedef std::vector<unsigned int>::const_iterator KeyIt;

    friend  ltkey;
    std::vector<unsigned int> m_key;
  };
}


        


                      
#endif
