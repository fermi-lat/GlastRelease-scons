// $Header$

#ifndef XMLUTIL_DICTFIELDMAN_H
#define XMLUTIL_DICTFIELDMAN_H

#include <string>

// g++ has hash_map; Visual Studio doesn't
#ifdef __GNUG__
#  include <hash_map>
#  define Registry hash_map<std::string, DictField*>

// otherwise might not have hash_map available, so use map, which
// is overkill for our situation since we don't need sorting
#else        
#  include <Map>
namespace {
  struct ltstr {
    bool operator()(const std::string s1, const std::string s2) const
    { return s1.compare(s2) < 0; }
  };
}

#  define Registry  map<std::string, DictField*, ltstr>
#else
#endif

namespace xmlUtil {
  //! Keep track of the collection of fields which may be referenced
  //! in a particular IdDict
  class DictField;

  class DictFieldMan {
  public:
    DictFieldMan(int size = 0);
    ~DictFieldMan();

    void register(DictField *field);
    const DictField const *find(std::string name);

  private:
    // Mechanism for keeping track of the fields is a 
    // pair-associative container with key = name and 
    // value = pointer to the DictField object
    Registry *reg;

  };
}
