// $Header$

#ifndef VolumeIdentifier_h
#define VolumeIdentifier_h
#include <deque>

namespace idents{

/**
Add some additional functionality to a vector of unsigned ints
Inherit from a deque for now. Should be private, but could not get
the append method to compile (?)
*/
class VolumeIdentifier : public std::deque<unsigned int>{
public:

    VolumeIdentifier();

    //! prepend, in front, another id
    void prepend( const VolumeIdentifier& id);

    //! append another id
    void append( const VolumeIdentifier& id);

    //! append an int
    void append(unsigned int id);

    //! return a name made up with slash delimiters
    std::string name(const char* delimiter="/")const ;

    
    //! return a long word for sorting 
#ifdef WIN32
    typedef __int64 int64;
#else
    typedef  long long int64;
#endif
    operator int64()const{ return m_value;}
    
private:
    void update(); // update value
    mutable int64 m_value; // for sorting
};

}
#endif