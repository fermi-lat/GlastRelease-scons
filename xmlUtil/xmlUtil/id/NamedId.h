//$Header$
#ifndef XMLUTIL_NAMEDID_H
#define XMLUTIL_NAMEDID_H

#include <vector>
#include <string>
#include <iostream>
#include "xmlUtil/id/Identifier.h"

namespace xmlUtil {

  //! An \b Identifier is simply a finite list of non-neg. integers
  //  typedef std::vector<unsigned> Identifier;

  //! Another convenience typedef for a list of strings (the
  //! id field names which may be extracted from a \b NamedId )
  typedef std::vector<std::string * > NameSeq;

  class NamedId;

  std::ostream& operator<<(std::ostream& s, const NameSeq& seq);
  std::ostream& operator<<(std::ostream& s, const Identifier& seq);
  std::ostream& operator<<(std::ostream& s, const NamedId& nId);

  class NamedId {
  public:
    friend std::ostream& operator<<(std::ostream& s, const NamedId& nId);

    //! When length is known, best to allocate all at once at start
    NamedId(const int len = 0);

    //! Copy constructor (should also have assignment operator)
    NamedId(const NamedId& toCopy);

    //! Just copy first n fields to a new NamedId
    NamedId(const NamedId& toCopy, unsigned int n);

    ~NamedId();

    struct IdField {
      std::string name;
      unsigned value;
    };

    //! Build a NamedId incrementally
    void addField(std::string name, unsigned value);

    void addField(const IdField& newField);

    unsigned int size() const {return m_fields->size();}
    //! Delete field on the end. 
    // Should it be returned rather than void?
    void popField(unsigned int n = 1);

    //! Return true if initial field names are compatible with \a subpath
    bool hasSubpath(const NameSeq& subpath) const;

    //! Return index of first occurrence of specified field, or -1
    //! if no occurrences.
    int hasField(std::string name) const;

    Identifier *stripNames();

    // Retrieve the sequence of names
    // So far this is just here for symmetry.  Can't think of a
    // good reason why anyone would want to do it.
    // NameSeq *stripValues() const;
    /* Return subsequence from field with name \b name1 to field with
      name \b name2, inclusive.  If these fields don't exist or don't
      occur in the proper order, return null.    */
    // NamedId *extract(const std::string& name1, const std::string& name2);

    /* Return a subsequence from field with name \b name1 until the
      end of the NamedId.  If \b name1 doesn't occur return null. */
    // NamedId *extract(const std::string& name1);

    // Append fields in order from \a toAppend 
    //void append(const NamedId& toAppend);

    friend class DictNode;
    friend class IdOperation;
    friend class IdOpCompress;

  protected:
    typedef std::vector <IdField* > Fields;
    typedef Fields::iterator  FieldIt;
   
    Fields *m_fields;
  };

}
#endif

