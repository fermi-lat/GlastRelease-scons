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
    packId(tower, vId[fLayer], vId[fCALLog], FACE_UNUSED, RANGE_UNUSED);
}

/*
  Meaning of vol id fields is  for CAL (see 
http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml#cal)
 is
    fLATObjects (must be 0 = eLATTowers)
    fTowerY     (range 0-3)
    fTowerX     (range 0-3)
    fTowerObjects (must be 0 = eTowerCAL)
    fLayer      (range 0-7)
    fMeasure      (0 for measures X, else 1)
    fCALLog     (range 0-11)
    fCellCmp    (0 for crystal, 1-4 for diodes)
    fCALSeg     (for crystal if segments are defined in geometry)

*/
VolumeIdentifier* CalXtalId::makeVolumeId() const {
  VolumeIdentifier* vid = new VolumeIdentifier;
  vid->append(0);   // LATObjects
  int fld = getTower();
  int towerY = fld/4;
  int towerX = fld - (4*towerY);
  vid->append(towerY);
  vid->append(towerX);
  vid->append(0);       // towerCAL
  fld = getLayer();
  vid->append(fld);
  if (isX()) {
    vid->append(0);
  } else {
    vid->append(1);
  }
  fld = getColumn();
  vid->append(fld);
  vid->append(0);   // only return vid's for crystals
  return vid;
}

// the inserter to stream unpacked tower, layer and column
void CalXtalId::write(std::ostream &stream) const
{
    stream << getTower() << " ";
    stream << getLayer() << " ";
    stream << getColumn() << " ";
    if (validFace() ) stream << getFace() << " ";
    if (validRange() ) stream << getRange() << " ";
}


// extract unpacked ID, and stuff packed ID from unpacked info
void CalXtalId::read(std::istream &stream)
{
    int tower,layer,column;
    stream >> tower >> layer >> column;
    packId(tower, layer, column, FACE_UNUSED, RANGE_UNUSED);		
}


// retrieve unpacked ID, tower, layer, and column
void CalXtalId::getUnpackedId(short& tower, short& layer, short& column)
{
    tower  = getTower();
    layer  = getLayer();
    column = getColumn();
}


