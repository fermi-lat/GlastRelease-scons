// $Header $

#include "xmlUtil/id/IdDictMan.h"
#include "xmlUtil/id/IdDict.h"
#include <algorithm>

namespace xmlUtil {
  IdDictMan * IdDictMan::me = 0;

  IdDictMan *IdDictMan::getPointer() {
    if (me == 0) {
      me = new IdDictMan();
    }
    return me;
  }

  IdDict *IdDictMan::findDict(std::string name) {
    unsigned int ix;

    for (ix = 0; ix < dictDir.size(); ix++) {
      IdDict *dict = (dictDir)[ix];
      if (name.compare(dict->getDictName()) == 0) return dict;
    }
    return 0;
  }
  
  IdDictMan::RetCode IdDictMan::registerDict(IdDict *dict) {
    // If already have a dictionary of this name, don't register
    if (findDict(dict->getDictName()) != 0) return DUPLICATE;

    dictDir.push_back(dict);
    return SUCCESS;
  }

  // Included for completeness.  
  IdDictMan::RetCode IdDictMan::removeDict(std::string name) {

    IdDict *dict = findDict(name);
    if (dict == 0) return NOTFOUND;

    std::remove(dictDir.begin(), dictDir.end(), dict);
    return SUCCESS;
  }


}
