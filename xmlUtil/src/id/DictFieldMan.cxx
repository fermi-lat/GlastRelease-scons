// $Header

#include "xmlUtil/id/DictFieldMan.h"
#include "xmlUtil/id/DictField.h"
#include <utility> 
#include <algorithm>


namespace xmlUtil {
  typedef pair<const char *, DictField*> RegPair;
  typedef Registry::iterator RegIterator;

  DictFieldMan::DictFieldMan(int size) {
    if (size > 0)   {
      reg = new Registry(size);
    }
    else   {
      reg = new Registry();
    }
    return;
  }

  DictFieldMan::~DictFieldMan() {
    // All other objects containing references to the fields 
    // (e.g., DictNodes) should have been deleted already.

    // Get rid of all fields...
    for (RegIterator it = reg->begin(); it != reg->end(); ++it) {
      
      delete it->second;
    }
    reg->clear();
    
      // then delete registry
    delete reg;
  }
  
  void DictFieldMan::signup(DictField *field)  {
    reg->insert(RegPair((field->getName()).c_str(), field));
  }

  const DictField* const DictFieldMan::find(const std::string& name) const {
    // find returns an iterator over the hash map, so a pointer
    // to a pair<std::string, DictField *>
    Registry::iterator  found;
    found = reg->find(name.c_str());
    if (found != reg->end() ) {
      return (*found).second;
    }
    else return 0;
      
      // *(reg->find(name));
  }

}   // end namespace
