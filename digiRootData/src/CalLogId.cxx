#include "digiRootData/CalLogId.h"

ClassImp(CalLogId)

// Include files
#include <iostream>
using namespace std;

#include "TObject.h"


/*!
//------------------------------------------------------------------------------
//
// ClassName:   CalLogId        
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
  //------------------------------------------------------------------------------
*/
CalLogId::CalLogId(Int_t packedId)
: m_packedId(packedId)
{
    unpackId();
}

CalLogId::CalLogId(Short_t tower, Short_t layer, Short_t column) 
{
    packId(tower, layer, column);
    m_tower = tower;
    m_layer = layer;
    m_column = column;
}

CalLogId::~CalLogId() {}

void CalLogId::init(){ }

// get packed ID or unpacked tower, layer, and column
Int_t CalLogId::getPackedId() const {return m_packedId;}
//void CalLogId::getUnpackedId(Short_t& tower, Short_t& layer, Short_t& column)

// get tower, layer, and column numbers individually from packed ID
Short_t CalLogId::getTower() const {return (m_packedId >> 0x7) & 0xf;}
Short_t CalLogId::getLayer() const {return (m_packedId >> 0x4) & 0x7;}
Short_t CalLogId::getColumn() const {return m_packedId & 0xf;}


/// retrieve unpacked ID, tower, layer, and column
void CalLogId::getUnpackedId(Short_t& tower, Short_t& layer, Short_t& column)
{
    tower  = m_tower;
    layer  = m_layer;
    column = m_column;
}

