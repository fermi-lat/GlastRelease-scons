// $Header$
#ifndef XMLUTIL_IDDICT_H
#define XMLUTIL_IDDICT_H

#include "dom/DOM_Element.hpp"

#include <vector>
#include <string>

#include "xmlUtil/id/NamedId.h"

namespace xmlUtil {
  //! IdDict corresponds to the xml element idDict in gdd.dtd
  class DictNode;
  class DictField;

  class IdDict {
  public:

    IdDict(DOM_Element elt);
    ~IdDict();

    //! Verify that dictionary is valid
    bool isValid() const;

    //! Verify that supplied identifier is valid w.r.t. the idDict.
    bool idOk(const Identifier& id) const;

    //! Verify that supplied named identifier is valid w.r.t. the idDict.
    bool idOk(const NamedId& id) const;

    //! Verify that the sequence of field names corresponds to at 
    //! least one valid identifier.
    bool nameSeqOk(const NameSeq& seq) const;

    //! Given an identifier (list of non-negative integers return
    //! the corresponding list of field names (empty list if invalid).
    // Not clear anyone wants or needs this service; if not, dump it
    //    NameSeq*     getNameSeq(const Identifier& id) const;

    //! Given an identifier, return a full NamedIdentifier (conceptually
    //! a sequence of (name, value) pairs).
    NamedId *  getNamedId(const Identifier& id) const;

    const std::string& getDictName() const {return m_name;};
    int                getMajorVersion() const {return m_major;};
    int                getMinorVersion() const {return m_minor;};
    int                getPatchVersion() const {return m_patch;};
    
  protected:
    friend class IdConverter;

    // Maybe make a derived class, IdDictBuildable, where the following
    // constructor and member function are declared.  IdConverter
    // would create an object of this (derived) type.

    //! Following constructor may be used by friend
    //! IdConverter to build output dictionary
    IdDict(DictNode* root);

    //! Link up newNode to parent; also update m_fields if necessary
    addChild(DictNode* parent, DictNode* newNode);

    vector<DictField*> m_fields; /*< all fields defined in the dictionary;
                                need not all be used (associated with
                                a node) */
    DictNode*          m_root;
    std::string        m_name;   /*< dictionary name */
    int                m_major;  /*< version information */
    int                m_minor;  /*< version information */
    int                m_patch;  /*< version information */
  };  // end of IdDict class definition
}      // end of namespace xmlUtil
#endif
