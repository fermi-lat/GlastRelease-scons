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

CalXtalId::CalXtalId(Short_t tower, Short_t layer, Short_t column) 
{
    packId(tower, layer, column);
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

void CalXtalId::init(Short_t tower, Short_t layer, Short_t column)
{ 
    Clear();
    packId(tower, layer, column);
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

Short_t CalXtalId::getColumn() const {
    return m_packedId & 0xf;
}

void CalXtalId::getUnpackedId(Short_t& tower, Short_t& layer, Short_t& column)
{
    tower  = getTower();
    layer  = getLayer();
    column = getColumn();
}

