// $Header$
#ifndef XMLUTIL_IDDICT_H
#define XMLUTIL_IDDICT_H

#include <dom/DOM_Element.hpp>

#include <vector>
#include <string>

#include "xmlUtil/id/DictObject.h"
#include "xmlUtil/id/NamedId.h"

namespace xmlUtil {
  //! IdDict corresponds to the xml element idDict in gdd.dtd
  /*! The structure of an IdDict is
       a collection of DictFields, each of which may point to an associated
       DictConstraints object.  The fields may be accessed via the single
       DictFieldMan object associated with the dictionary
       A collection of nodes arranged hierarchically.  Each node
          MUST reference its parent node (unless it's the root
            of the hierarchy, in which case it has none)
          MUST reference a DictField object in the collection
          MAY reference an additional DictConstraints object which
            further restricts the values its associated field may
            take on
          MAY referemce yet another DictConstraints object which
            restricts the values its parent node field may take on
  */

  /* Some day should redo this using Builder architecture.
     Define abstract (or nearly abstract) IdDictBuilder class, then
     concrete derived XMLbuilder.  That way the dictionary class and
     classes making up elements of the dictionary wouldn't have to know
     about the XML representation at all.
  */
  class DictNode;
  class DictField;
  class DictFieldMan;
  class DictVisitor;

  class IdDict : public DictObject {
  public:
    enum Constituents {nodeHierarchy = 1,
                       fieldManager = 2};

    IdDict(DOM_Element elt);
    ~IdDict();

    //! Verify that dictionary is valid
    bool isValid();

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
    DictNode&          getRoot() const {return *m_root;};
    bool               accept(DictVisitor *vis);

    bool               accept(DictVisitor *vis, unsigned constituentMask);
    //! Extra accept function to allow caller to choose among
    //! constituents of dictionary to visit.  May choose all,
    //! just node hierarchy, or just field manager.
  protected:
    friend class IdConverter;

    // Maybe make a derived class, IdDictBuildable, where the following
    // constructor and member function are declared.  IdConverter
    // would create an object of this (derived) type.

    //! Following constructor may be used by friend
    //! IdConverter to build output dictionary
    IdDict(DictNode* root);

    //! Link up newNode to parent; also ask m_fieldMan to add the
    //! field if it's new
    bool addChild(DictNode* parent, DictNode* newNode);

    /* vector<DictField*> m_fields;  */
                            /*< all fields defined in the dictionary;
                                need not all be used (associated with
                                a node) */
    DictFieldMan*      m_fieldMan;
    DictNode*          m_root;
    std::string        m_name;   /*< dictionary name */
    int                m_major;  /*< version information */
    int                m_minor;  /*< version information */
    int                m_patch;  /*< version information */
  };  // end of IdDict class definition
}      // end of namespace xmlUtil
#endif
