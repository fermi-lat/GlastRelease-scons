// $Header$
#ifndef XMLUTIL_NAMEDID_H
#define XMLUTIL_NAMEDID_H

#include <vector>
#include <string>

namespace xmlUtil {

  //! An /b Identifier is simply a finite list of non-negative
  //! integers.  
  typedef std::vector<unsigned> Identifier;

  //! Another convenience typedef for a list of strings (the
  //! id field names which may be extracted from a /b NamedId )
  typedef std::vector<std::string*> NameSeq;

  class NamedId {
  public:

    //! Default constructor creates empty NamedId
    //  (Perhaps allocate a vector of modest size to start??)
    NamedId() {};

    //! More efficient to allocate resources all at once when
    //! final length of NamedId is known in advance.
    NamedId(const int len);

    //! Copy constructor
    NamedId(const NamedId& toCopy);

    ~NamedId() {};


    //! Build a NamedId incrementally
    void addField(const std::string name, const unsigned value);
    
    //! Retrieve just the value of each field
    Identifier *extractIdentifier() const;

    //! Retrieve the sequence of names
    // So far this is just here for symmetry -- can't think of a 
    // good reason why anyone would want to do it.
    NameSeq    *extractNameSeq() const;

  private:
    typedef struct s_IdField {
      std::string name;
      unsigned value;
    }  IdField;

    typedef vector<*IdField> Fields;

    Fields *m_fields;
  };
}
#endif
