
#include "mcRootData/VolumeIdentifier.h"
#include <iostream>

ClassImp(VolumeIdentifier)

// File and Version Information:
//      $Header$
//
// Description:
//      The class VolumeIdentifier encapsulates volume identifiers defined in
//      the xml file describing the detector geometry. It represents the 
//      identifier through two 32 integers to enable an efficient sorting of 
//      objects which use this identifier.
//      Every single identifier which constitute the volume identifier is code
//      into a binary string of 6 bits and these bit strings are packed into two
//      32 bit integers.  So, every single id can be an integer beteween 0 and 
//      63 and each volume identifier can be built by a maximum number of 10 ids.
//      index == 0 corresponds to the upper most bits  
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

VolumeIdentifier::VolumeIdentifier(const VolumeIdentifier& id) {
    m_bits0to31 = id.m_bits0to31;
    m_bits32to63 = id.m_bits32to63;
    m_size = id.m_size;
}

void VolumeIdentifier::initialize(UInt_t bits0to31, UInt_t bits32to63, UInt_t size)
{
    m_bits0to31 = bits0to31;
    m_bits32to63 = bits32to63;
    m_size = size;
}

void VolumeIdentifier::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "Size: " << m_size << endl;
    cout << this->name() << endl;
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
    
    unsigned int i;
    s << delimiter;
    for (i = 0; i < m_size; i++)
    {
        s << this->operator[](i) << delimiter;
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

UInt_t VolumeIdentifier::operator[](unsigned int index) const
{
    // Purpose and Method:  Return a particular id corresponding to the index
    //    Note that index == 0 refers to the upper most bits, in m_bits32to63

    // Check for invalid index, we only store 10 ids per VolumeIdentifier
    if (index > 9) return 0;

    // Turn on 6 bits:  111111
    static UInt_t mask_id = 63;
    const UInt_t numBitsPerId = 6;
    const UInt_t shift = 24;

    // Stored the first 5 ids (index 0-4) in m_bits32to63
    if (index < 5) {
        UInt_t idBits = m_bits32to63 >> (shift - numBitsPerId*index);
        return (idBits & mask_id);
    } 
    // other 5 ids are in m_bits0to31
    UInt_t idBits = m_bits0to31 >> (shift - 6*(index-5));
    return (idBits & mask_id);
}


void VolumeIdentifier::append( unsigned int id)
{
    // Purpose and Method:  Add another 6 bit id to the packed bit string
    // the first id appended becomes the most significant number in the internal
    // representation. In this way I can obtain an equivalent of the lexicographic order
    // between volume identifiers
    const UInt_t numBitsPerId = 6;
    const UInt_t shift = 24;

    UInt_t t = id;
    // We store the first 5 ids (index 0-4) in m_bits32to63
    if (m_size < 5) {
        t = t << shift; 
        // Filling the upper bits
        m_bits32to63 = m_bits32to63 | (t >> (numBitsPerId * m_size));
    } else {
        // Store the second 5 ids (idex 5-9) in m_bits0to31
        t = t << shift;
        m_bits0to31 = m_bits0to31 | (t >> numBitsPerId * (m_size-5));
    } 
    
    m_size++;
}
