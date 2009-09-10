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

#include <sstream>
#include <cassert>
#include <stdexcept>

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
    std::stringstream s;
    
    unsigned int bufIds = 0;
    static unsigned int mask = (1 << s_bitsPer) - 1; //  6 bits for final mask
   
    int64 copyValue = m_value;
    
    int i;
    s << delimiter;

    // Fields are written out starting with one located in most significant
    // bits.  Since 64 is not evenly divisible by s_bitsPer (6) this field
    // is located in bits 54-59.  Top 4 bits are unused.
    for (i = 0; i < m_size; i++)
    {
        bufIds = ((copyValue ) >> s_maxShift) & mask;
        s << bufIds << delimiter;
        copyValue = copyValue << s_bitsPer;
    }
    
    std::string tmp=s.str();
    return tmp.substr(0,tmp.size()-1);
}


unsigned int VolumeIdentifier::operator[](unsigned int index)
{
    static int64 mask2 = (1 << s_bitsPer) - 1;            /* 63 */
    int64 copyShifted = m_value >> (s_maxShift - s_bitsPer*index);
    return (copyShifted & mask2);
}

unsigned int VolumeIdentifier::operator[](unsigned int index) const
{
  static int64 mask2 = (1 << s_bitsPer) -1 ; // 63;
    int64 copyShifted = m_value >> (s_maxShift - s_bitsPer*index);
    return (copyShifted & mask2);
}


void VolumeIdentifier::prepend( const VolumeIdentifier& id)
{
    m_value = (m_value >> (s_bitsPer*id.size())) | id.getValue();
    m_size += id.size();
}

void VolumeIdentifier::append( const VolumeIdentifier& id)
{
    m_value = m_value | (id.getValue() >> (s_bitsPer*m_size));
    m_size += id.size();
}

void VolumeIdentifier::append( unsigned int id)
{
  // the first id appended becomes the most significant digit in the internal
  // rappresentation. In this way I can obtain an equivalent of the lexicographic order
  // between volume identifiers
  if (m_size >= (int) s_maxSize) {
    std::string 
      errtxt("VolumeIdentifier::append: id is already of maximum size");
    throw std::range_error(errtxt);
  }
  else if (id > s_maxFieldValue) {
    std::string 
      errtxt("VolumeIdentifier::append: new field value is too large");
    throw std::range_error(errtxt);
  }
  int64 t = id;
  t &= (1 << s_bitsPer) - 1;    // only use s_bitsPer bits from id.
  t = t << s_maxShift;
  m_value = m_value | (t >> s_bitsPer*m_size);
  m_size++;
}
