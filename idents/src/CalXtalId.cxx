// File and Version information
// $Header$
//
// ClassName:   CalXtalId        
//  
// Description: ID class for CAL Xtals                                
//		Supports both packed Xtal ID = (tower*8 + layer)*16 + column
//		 and	unpacked ID, i.e. tower, layer, and column.
//		Extractor/inserter give i/o for unpacked ID.
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



// Include files
#include "idents/CalXtalId.h"
#include "idents/VolumeIdentifier.h"
#include <stdexcept>

using namespace idents; 

// Constructor from VolumeIdentifier.  Assuming a particular form for
// fields within VolumeIdentifier (correct as of June 14, 2004):
//    field 0 is fLATObjects; check for value = 0 (towers)
//    field 1 is y-value for tower 
//    field 2 is x-value for tower
//    feild 3 is fTowerObjects; check for value = 0 (Calorimeter)
//    field 4 is Cal layer
//    field 5 is orientation (measures X or Y)
//    field 6 log number ("column" in CalXtalId terms)
CalXtalId::CalXtalId(const VolumeIdentifier& vId, unsigned xNum) {
    const int minSize = 7;
    const unsigned LATObjectTower = 0;
    const unsigned TowerObjectCal = 0;
    enum {
        fLATObjects = 0,
        fTowerY = 1,
        fTowerX = 2,
        fTowerObjects = 3,
        fLayer = 4,
        fMeasure = 5, 
        fCALLog = 6
    };
    if (vId.size() < minSize) throw std::invalid_argument("VolumeIdentifier");
    if ((vId[fLATObjects] != LATObjectTower) || 
        (vId[fTowerObjects] != TowerObjectCal)) {
        throw std::invalid_argument("VolumeIdentifier");
    }
    if (vId[fTowerX] >= xNum) throw std::invalid_argument("xNum");
    short tower = xNum*vId[fTowerY] + vId[fTowerX];
    packId(tower, vId[fLayer], vId[fCALLog]);
}

// the inserter to stream unpacked tower, layer and column
void CalXtalId::write(std::ostream &stream) const
{
    stream << getTower() << " ";
    stream << getLayer() << " ";
    stream << getColumn() << " ";
}


// extract unpacked ID, and stuff packed ID from unpacked info
void CalXtalId::read(std::istream &stream)
{
    int tower,layer,column;
    stream >> tower >> layer >> column;
    packId(tower, layer, column);		
}


// retrieve unpacked ID, tower, layer, and column
void CalXtalId::getUnpackedId(short& tower, short& layer, short& column)
{
    tower  = getTower();
    layer  = getLayer();
    column = getColumn();
}


