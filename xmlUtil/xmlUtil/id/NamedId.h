//$Header$
#ifndef XMLUTIL_NAMEDID_H
#define XMLUTIL_NAMEDID_H

#include <vector>
#include <string>
namespace xmlUtil {

  //! An \b Identifier is simply a finite list of non-neg. integers
  typedef std::vector<unsigned> Identifier;

  //! Another convenience typedef for a list of strings (the
  //! id field names which may be extracted from a \b NamedId )
  typedef std::vector<std::string * > NameSeq;

  class NamedId {
  public:


    //! When length is known, best to allocate all at once at start
    NamedId(const int len = 0);

    //! Copy constructor (should also have assignment operator)
    NamedId(const NamedId& toCopy);

    ~NamedId();

    //! Build a NamedId incrementally
    void addField(const std::string name, const unsigned value);

    //! Return true if initial field names are compatible with \a path
    bool hasPath(const NamedSeq& path);

    bool hasField(const std::string name);

    // Retrieve just the value of each field
    // Identifier *extractIdentifier() const;

    // Retrieve the sequence of names
    // So far this is just here for symmetry.  Can't think of a
    // good reason why anyone would want to do it.
    // NameSeq *extractNameSeq() const;

  private:
    friend class DictNode;

    typedef struct s_IdField {
      std::string name;
      unsigned value;
    }  IdField;

    typedef std::vector <IdField* > Fields;

    typedef Fields::iterator  FieldIt;
    
    Fields *m_fields;
  };
}
#endif

