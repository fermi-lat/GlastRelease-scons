// $Header$

#ifndef XMLUTIL_DICTFIELDMAN_H
#define XMLUTIL_DICTFIELDMAN_H

#include <string>

// g++ has hash_map; Visual Studio doesn't
#ifdef __GNUG__
#  include <hash_map>
#  define Registry hash_map<const char *, DictField*>
// otherwise might not have hash_map available, so use map, which
// is overkill for our situation since we don't need sorting
#else        
#  include <Map>
#  define Registry  map<const char *, DictField*>
#endif

namespace xmlUtil {
  //! Keep track of the collection of fields which may be referenced
  //! in a particular IdDict
  class DictField;

  class DictFieldMan {
  public:
    DictFieldMan(int size = 0);
    ~DictFieldMan();

    void signup(DictField *field);
    const DictField * const find(const std::string& name) const;

  private:
    // Mechanism for keeping track of the fields is a 
    // pair-associative container with key = name and 
    // value = pointer to the DictField object
    Registry *reg;

  };
}
#endif
