#include "commonRootData/idents/CalXtalId.h"

ClassImp(CalXtalId)

#include <iostream>
using namespace std;

#include "TObject.h"

//
// ClassName:   CalXtalId        
//  
// Description: ID class for CAL logs based on CalXtalId in the idents package                               
//		Supports both packed log ID = (tower*8 + layer)*16 + column and
//			unpacked ID, i.e. tower, layer, and column.
//
//              
// Author:  Heather Kelly - based on idents::CalXtalId written by J. Eric Grove

CalXtalId::CalXtalId() {
    Clear();
}

CalXtalId::CalXtalId(UInt_t packedId)
: m_packedId(packedId)
{
}

CalXtalId::CalXtalId(Short_t tower, Short_t layer, Short_t column,
                     Short_t face, Short_t range) 
{
    packId(tower, layer, column, face, range);
}

CalXtalId::~CalXtalId() {

}

void CalXtalId::Clear(Option_t *option) {
    m_packedId = 0;
}

void CalXtalId::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Packed Id: " << m_packedId << endl;
    cout << "(Tower, Layer, Col): " << "(" << getTower() 
        << "," << getLayer() << "," << getColumn() << ")" << endl;
}

void CalXtalId::init(Short_t tower, Short_t layer, Short_t column, 
                     Short_t face, Short_t range)
{ 
    Clear();
    packId(tower, layer, column, face, range);
}

void CalXtalId::init(UInt_t packedId) {
    m_packedId = packedId;
}

UInt_t CalXtalId::getPackedId() const { 
    return m_packedId;
}

Short_t CalXtalId::getTower() const {
    return (m_packedId >> 0x7) & 0xf;
}

Short_t CalXtalId::getLayer() const {
    return (m_packedId >> 0x4) & 0x7;
}

Short_t CalXtalId::getFace() const {
    if (!validFace()) return FACE_UNUSED;
    return ((m_packedId >> FACE_SHIFT) & 0x1);
}
Short_t CalXtalId::getRange() const {
    if (!validRange()) return RANGE_UNUSED;
    return ((m_packedId >> RANGE_SHIFT) & 0x3);
}

Short_t CalXtalId::getColumn() const {
    return m_packedId & 0xf;
}

void CalXtalId::getUnpackedId(Short_t& tower, Short_t& layer, Short_t& column)
{
    tower  = getTower();
    layer  = getLayer();
    column = getColumn();
}

