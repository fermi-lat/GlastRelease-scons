#include "digiRootData/CalXtalId.h"

ClassImp(CalXtalId)

#include <iostream>
using namespace std;

#include "TObject.h"

//
// ClassName:   CalXtalId        
//  
// Description: ID class for CAL logs                                
//		Supports both packed log ID = (tower*8 + layer)*16 + column and
//			unpacked ID, i.e. tower, layer, and column.
//		Serialize methods give i/o for packed ID.
//		Extractor/inserter give i/o for unpacked ID.
//		Input methods ensure that packed and unpacked IDs stay in sync.
//
//    Retrieve packed ID or unpacked tower, layer, and column 
//        inline int getPackedId() const 
//        void getUnpackedId(short& tower, short& layer, short& column); 
// 
//    Retrieve tower, layer, and column numbers individually from packed ID 
//        inline short getTower() const 
//        inline short getLayer() const 
//        inline short getColumn() const 
//              
// Author:  J. Eric Grove	22 Mar 2001                      


CalXtalId::CalXtalId(UInt_t packedId)
: m_packedId(packedId)
{
    unpackId();
}

CalXtalId::CalXtalId(Short_t tower, Short_t layer, Short_t column) 
{
    packId(tower, layer, column);
    m_tower = tower;
    m_layer = layer;
    m_column = column;
}

CalXtalId::~CalXtalId() {

}

void CalXtalId::Clear(Option_t *option) {
    m_tower = 0;
    m_layer = 0;
    m_column = 0;
    m_packedId = 0;
}

void CalXtalId::Print(Option_t *option) const {

}

void CalXtalId::init(Short_t tower, Short_t layer, Short_t column)
{ 
    Clear();
    m_tower = tower;
    m_layer = layer;
    m_column = column;
    packId(tower, layer, column);
}

void CalXtalId::init(UInt_t packedId) {
    m_packedId = packedId;
    unpackId();
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
    tower  = m_tower;
    layer  = m_layer;
    column = m_column;
}

void CalXtalId::unpackId() 
{
    m_tower = getTower();
    m_layer = getLayer();
    m_column = getColumn();
}

