// File and Version Information:
//      \$Header\$
//
// Description:
//      The class VolumeIdentifier encapsulates volume identifiers defined in the xml
//      file describing the detector geometry. It rappresents the identifier through
//      a 64 integer to enable an efficient sorting of objects which use this identifier.
//      Every single identifier which constitutes the volume identifier is coded into a 
//      binary string of 6 bits and these bit strings are packed into a 64 bit integer.
//      So, every single id can be an integer beteween 0 and 63 and each volume identifier
//      can be built by a maximum number of 10 ids.  
//
// Author(s):
//      Toby Burnett		
//      Marco Frailis		


#include "idents/VolumeIdentifier.h"

#ifndef WIN32  
#include <strstream>
#else
#include <sstream>
#endif
#include <cassert>

using namespace idents;


VolumeIdentifier::VolumeIdentifier():  m_value(0), m_size(0){}

void VolumeIdentifier::init(int64 value, unsigned int size)
{
    m_value = value;
    m_size = size;
}

// Return the equivalent string of the volume identifier, that is the single
// ids separated by a '/' character
std::string VolumeIdentifier::name(const char* delimiter) const
{
#ifndef WIN32    
    std::strstream s;
#else
    std::stringstream s;
#endif
    
    unsigned int bufIds = 0;
    
    // this is a 64 bit mask with 6 bits (positions 54-59) set to 1 and the others to 0
    static int64 mask = 0x0fc0000000000000;
    
    int64 copyValue = m_value;
    
    int i;
    s << delimiter;
    for (i = 0; i < m_size; i++)
    {
        bufIds = (copyValue & mask) >> 54;
        s << bufIds << delimiter;
        copyValue = copyValue << 6;
    }
    
#ifndef WIN32
    s << '\0';
#endif
    std::string tmp=s.str();
    return tmp.substr(0,tmp.size()-1);
}


unsigned int VolumeIdentifier::operator[](unsigned int index)
{
    static int64 mask2 = 63;
    int64 copyShifted = m_value >> (54 - 6*index);
    return (copyShifted & mask2);
}

void VolumeIdentifier::prepend( const VolumeIdentifier& id)
{
    m_value = (m_value >> (6*id.size())) | id.getValue();
    m_size += id.size();
}

void VolumeIdentifier::append( const VolumeIdentifier& id)
{
    m_value = m_value | (id.getValue() >> (6*m_size));
    m_size += id.size();
}

void VolumeIdentifier::append( unsigned int id)
{
    // the first id appended becomes the most significant digit in the internal
    // rappresentation. In this way I can obtain an equivalent of the lexicographic order
    // between volume identifiers
    int64 t = id;
    t = t << 54;
    m_value = m_value | (t >> 6*m_size);
    m_size++;
}
