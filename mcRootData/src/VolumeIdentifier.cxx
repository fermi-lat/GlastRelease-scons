
#include "mcRootData/VolumeIdentifier.h"

ClassImp(VolumeIdentifier)

// File and Version Information:
//      $Header$
//
// Description:
//      The class VolumeIdentifier encapsulates volume identifiers defined in
//      the xml file describing the detector geometry. It represents the 
//      identifier through a 64 integer to enable an efficient sorting of 
//      objects which use this identifier.
//      Every single identifier which constitute the volume identifier is code
//      into a binary string of 6 bits and these bit strings are packed into a
//      64 bit integer.  So, every single id can be an integer beteween 0 and 
//      63 and each volume identifier can be built by a maximum number of 10 ids.  
//


#ifndef WIN32  
#include <strstream>
#else
#include <sstream>
#endif


VolumeIdentifier::VolumeIdentifier() :  
m_bits0to31(0), m_bits32to63(0), m_size(0)
{
}

void VolumeIdentifier::initialize(UInt_t bits0to31, UInt_t bits32to63, UInt_t size)
{
    m_bits0to31 = bits0to31;
    m_bits32to63 = bits32to63;
    m_size = size;
}

std::string VolumeIdentifier::name(const char* delimiter) const
{
    // Purpose and Method:  Return the equivalent string of the volume 
    //    identifier, that is the single ids separated by a '/' character
#ifndef WIN32    
    std::strstream s;
#else
    std::stringstream s;
#endif
    
    UInt_t bufIds = 0;
    
    static UInt_t mask_bits0to31 = 0xfc000000;
    // this is a 64 bit mask with the 6 bits (positions 54-59) set to 1 and the
    // others to 0
    static UInt_t mask_bits32to63 = 0x0fc00000;
    
    UInt_t copy_bits0to31 = m_bits0to31;
    UInt_t copy_bits32to63 = m_bits32to63;
    
    unsigned int i;
    s << delimiter;
    for (i = 0; i < m_size; i++)
    {
        if (i < 4) {
            bufIds = (copy_bits32to63 & mask_bits32to63) >> 22;
            s << bufIds << delimiter;
            copy_bits32to63 = copy_bits32to63<< 6;
        } else if (i == 4) {
            UInt_t mask_upper4bits = 15 << 24;
            UInt_t mask_lower2bits = 3 << 30;
            UInt_t upperBits = (copy_bits32to63 && mask_upper4bits) >> 22;
            UInt_t lowerBits = (copy_bits0to31 && mask_lower2bits) >> 30;
            bufIds = upperBits | lowerBits;
            s << bufIds << delimiter;
            copy_bits0to31 = copy_bits0to31 << 2;
            
        } else {
            bufIds = (copy_bits0to31 & mask_bits0to31) >> 26;
            s << bufIds << delimiter;
            copy_bits0to31 = copy_bits0to31 << 6;  
        }
    }
    
#ifndef WIN32
    s << '\0';
#endif
    std::string tmp=s.str();
    return tmp.substr(0,tmp.size()-1);
}

void VolumeIdentifier::Clear(Option_t *option) {
    m_bits0to31 = 0;
    m_bits32to63 = 0;
    m_size = 0;
}

UInt_t VolumeIdentifier::operator[](unsigned int index)
{
    // Check for invalid index, we only store 10 ids per VolumeIdentifier
    if (index > 9) return 0;
    static UInt_t mask_id = 63;
    
    if (index < 4) {
        UInt_t idBits = m_bits32to63 >> (22 - 6*index);
        return (idBits & mask_id);
    } else if (index == 4) {
        UInt_t mask_lower4bits = 15;
        UInt_t mask_upper2bits = 3 << 30;
        UInt_t upperBits = m_bits32to63 && mask_lower4bits;
        UInt_t lowerBits = m_bits0to31 && mask_upper2bits;
        UInt_t id = (upperBits << 2) | lowerBits;
        return id;
    } 
    UInt_t idBits = m_bits0to31 >> (22 - 6*(index-5));
    return (idBits & mask_id);
}


void VolumeIdentifier::append( unsigned int id)
{
    // the first id appended becomes the most significant number in the internal
    // representation. In this way I can obtain an equivalent of the lexicographic order
    // between volume identifiers
    UInt_t t = id;
    if (m_size < 4) {
        t = t << 22; //54;
        // Filling the upper bits
        m_bits32to63 = m_bits32to63 | (t >> 6 * m_size);
    } else if (m_size == 4) {
        // Splitting bits across the two 32 bit ints, 4 bits in the upper, 2 bits in the lower
        UInt_t copyT = t << 30;
        t = t << 22; //54;
        m_bits32to63  = m_bits32to63 | (t >> 6 * m_size);
        m_bits0to31 = m_bits0to31 | copyT;
    } else {
        t = t << 24;
        // Storing in the lower bits
        m_bits0to31 = m_bits0to31 | (t >> 6*(m_size-5));
    } 
    
    m_size++;
}
