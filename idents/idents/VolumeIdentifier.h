
#ifndef VolumeIdentifier_h
#define VolumeIdentifier_h

#include <string>

/** 
 * @class VolumeIdentifier
 *
 * @brief This class encapsulates volume identifiers.
 *
 * This class encapsulates volume identifiers defined in the xml files describing 
 * the detector. The volume identifiers are transformed into a 64 bit integer for 
 * efficiency purpose.
 *  
 * @author Toby Burnett
 * @author Marco Frailis
 *    
 * \$Header\$
 */

namespace idents{


class VolumeIdentifier {
public:

#ifdef WIN32
    typedef __int64 int64;
#else
    typedef  long long int64;
#endif

    VolumeIdentifier();

    /** 
     * This method initialize the VolumeIdentifier with a 64 bit integer and a size
     * and is used by the overloaded operator<< to build an identifier reading from
     * a persistent data store
     */
    init(int64 , unsigned int); 

    /// prepend, in front, another id
    void prepend( const VolumeIdentifier& id);

    /// append another id
    void append( const VolumeIdentifier& id);

    /// append an int
    void append(unsigned int id);

    /// number of single ids which constitute the volume identifier
    int size() const { return m_size;}

    /// return a name made up with slash delimiters
    std::string name(const char* delimiter="/")const ;

    
    /// return a long word for sorting 

    operator int64()const{ return m_value;}

    /// access single ids which constitute the volume identifier
    unsigned int operator[](unsigned int);
    
private:

    /// internal rappresentation of the volume identifier
    mutable int64 m_value; // for sorting
    /// number of single ids which constitute the volume identifier
    int m_size;
};

}
#endif
