//$Header$
#ifndef XMLUTIL_IDENTIFIER_H
#define XMLUTIL_IDENTIFIER_H

#include <vector>

namespace xmlUtil {
  /// Class used primarily to incrementally build identifiers;
  /// see also xmlUtil::IdKey.
  /** For now, nothing to add to base class. */
  class Identifier : public std::vector<unsigned int> {
  public:
    Identifier();
    ~Identifier() {}
    Identifier(unsigned n);
    Identifier(const Identifier& toCopy);
    void append(unsigned int n) { push_back(n);}
    void append(const Identifier& toAppend);
  };
}
#endif
