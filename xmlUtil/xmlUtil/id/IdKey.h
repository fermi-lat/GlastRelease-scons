// $Header$
#ifndef XMLUTIL_IDKEY_H
#define XMLUTIL_IDKEY_H

#include <vector>
#include <deque>
#include <string>
#include <iterator>


namespace xmlUtil {

  ///  IdKey is a compressed, sortable version of an identifier which may
  ///  be used as a key for a map.
  class IdKey {
  public:
    /// Construct an IdKey from  a vector
    IdKey(const std::vector<unsigned int> id);

    /// Construct an IdKey from  a deque
    IdKey(const std::deque<unsigned int> id);

    IdKey() {}

    /*
       Construct an IdKey from an identifier whose fields are expressed
       as strings         NYI
          IdKey(const std::vector<std::string>&);  
       Do we also want to go the other direction? (IdKey to std::string)
   */

    ~IdKey();

    /// The comparison function needed when making a map
    class ltkey { 
    public:
      bool operator()(IdKey key1, IdKey key2) const;
    };
    friend  class ltkey;       

  private:
    /** template function which does the actual work of building
        the new IdKey.  In principle it could be instantiated for
        any container whose value type is unsigned int. */
    template<class Container>
    void fillFrom(const Container& con);

    typedef std::vector<unsigned int>::const_iterator KeyIt;
    /// Compressed identifier is stored in /a m_key
    std::vector<unsigned int> m_key; 
  };
}


        


                      
#endif
