// $Header$

#include "xmlUtil/id/DictValidVisitor.h"
// #include "xmlUtil/IdDict.h"
#include "xmlUtil/id/DictNode.h"

namespace xmlUtil {

  bool DictValidVisitor::visitDict() {
    if (!m_dict) {
      m_state = eInvalid;
      return false;
    }
    return visitDict(m_dict);
  }

  bool DictValidVisitor::visitDict(IdDict* dict) {
    m_dict = dict;

    m_state = eValid;  // so far

    return true;
  }

  bool DictValidVisitor::visitNode(DictNode* node) {
    if (m_state != eValid) {
      m_state = eInvalid;
      return false;
    }
    if ((node->consistentChildren()) &&
        (node->consistentParent())   &&
        (node->consistentValues())     ) return true;

    m_state = eInvalid;
    return false;
  }

}  // end namespace
