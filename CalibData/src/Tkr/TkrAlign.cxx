// $Header$

/// Implementations for all tracker-alignment related calib. TDS classes

#include "CalibData/Tkr/TkrTowerAlignCalib.h"
#include "CalibData/Tkr/TkrInternalAlignCalib.h"
#include "CalibData/CalibModel.h"
#include <stdexcept>

using CLHEP::Hep3Vector;

class MsgStream;

// Make a key for tower, tray, face, ladder or wafer as follows:
// Set lower 4 bits to indicate which fields are valid.  Tower is
// always valid
// Next 2 bits (4-5) represent wafer
// 6-7 represent ladder
// 8 represents face
// 9-13 are for tray
// 14-17 are for tower
#define HAS_WAFER 1
#define HAS_LADDER 2
#define HAS_FACE  4
#define HAS_TRAY  8
#define WAFER_SHIFT 4
#define LADDER_SHIFT 6
#define FACE_SHIFT 8
#define TRAY_SHIFT 9
#define TOWER_SHIFT 14


#define TKR_MAX_TOWER_ID 15
#define TKR_MAX_TRAY_ID  31
#define TKR_MAX_FACE_ID 1
#define TKR_MAX_LADDER_ID 3
#define TKR_MAX_WAFER_ID 3


namespace {
  void assign3(const Hep3Vector& inVec, Hep3Vector& outVec) {
    for (unsigned ix = 0; ix < 3 ; ix++) outVec[ix] = inVec[ix];
  }
}

namespace CalibData {
  
  void TkrTowerAlign::getAlign(CLHEP::Hep3Vector& disp, CLHEP::Hep3Vector& rot)
  {
    assign3(m_disp, disp);
    assign3(m_rot, rot);
  }

  void TkrTowerAlign::putAlign(const CLHEP::Hep3Vector& disp, 
                               const CLHEP::Hep3Vector& rot) {
    assign3(disp, m_disp);
    assign3(rot, m_rot);
  }

  void TkrTowerAlign::update(RangeBase* other) {
    // if dynamic cast doesn't work something is very wrong
    TkrTowerAlign* other1 = dynamic_cast<TkrTowerAlign*> (other);
    std::string msg("TkrTowerAlign::update");
    if (!other1) throw std::invalid_argument(msg);
    assign3(other1->m_rot, m_rot);
    assign3(other1->m_disp, m_disp);
  }

  StatusCode TkrTowerAlignCalib::update(CalibBase& other, 
                                        MsgStream*  log ) {
    TkrTowerAlignCalib& other1 = dynamic_cast<TkrTowerAlignCalib& > (other);
    StatusCode sc = CalibBase::update(other1, log);
    if (sc != StatusCode::SUCCESS) return sc;

    unsigned sz = other1.m_towers.size();

    if (sz != m_towers.size()) return StatusCode::FAILURE;
    for (unsigned ix = 0; ix < sz; ix++) {
      Hep3Vector disp, rot;
      (other1.m_towers[ix]).getAlign(disp, rot);
      m_towers[ix].putAlign(disp, rot);
    }
    return StatusCode::SUCCESS;
  }


  TkrInternalAlign::TkrInternalAlign() {
    assign3(Hep3Vector(0,0,0), m_disp);
    assign3(Hep3Vector(0,0,0), m_rot);
  }


  TkrInternalAlignCalib::~TkrInternalAlignCalib() {
    clean();
  }

  void TkrInternalAlignCalib::clean() {
    AlignMap::iterator it = m_items.begin();
    while (it != m_items.end() ) {
      if (it->second != 0) delete it->second;
      ++it;
    }
    m_items.clear();
  }
    

  const CLID&  TkrTowerAlignCalib::classID() {
    return CLID_Calib_TKR_TowerAlign;
  }


  const CLID&  TkrInternalAlignCalib::classID() {
    return CLID_Calib_TKR_InternalAlign;
  }

  void TkrInternalAlign::getAlign(CLHEP::Hep3Vector& disp, 
                                  CLHEP::Hep3Vector& rot)  {
    assign3(m_disp, disp);
    assign3(m_rot, rot);
  }

  void TkrInternalAlign::putAlign(const CLHEP::Hep3Vector& disp, 
                                  const CLHEP::Hep3Vector& rot) {
    assign3(disp, m_disp);
    assign3(rot, m_rot);
  }



  unsigned TkrInternalAlignCalib::makeTowerKey(unsigned tower) const {
    return tower << TOWER_SHIFT;
  }
  unsigned TkrInternalAlignCalib::makeTrayKey(unsigned tower, unsigned tray) 
  const {
    return (tower << TOWER_SHIFT) + (tray << TRAY_SHIFT) + HAS_TRAY;
  }
  unsigned TkrInternalAlignCalib::makeFaceKey(unsigned tower, unsigned tray, 
                                              unsigned face) const {
    return (tower << TOWER_SHIFT) + (tray << TRAY_SHIFT) + 
      (face << FACE_SHIFT) + HAS_FACE;
  }
  unsigned TkrInternalAlignCalib::makeLadderKey(unsigned tower, unsigned tray, 
                                                unsigned face, unsigned ladder)
  const {
    return (tower << TOWER_SHIFT) + (tray << TRAY_SHIFT) + 
      (face << FACE_SHIFT) + (ladder << LADDER_SHIFT) + HAS_LADDER;
  }
  unsigned TkrInternalAlignCalib::makeWaferKey(unsigned tower, unsigned tray, 
                                               unsigned face, unsigned ladder, 
                                               unsigned wafer) const {
    return (tower << TOWER_SHIFT) + (tray << TRAY_SHIFT) + 
      (face << FACE_SHIFT) + (ladder << LADDER_SHIFT) +
      (wafer << WAFER_SHIFT) + HAS_WAFER;
  }


  bool
  TkrInternalAlignCalib::checkLimits(unsigned* tower, unsigned* tray, 
                                     unsigned* face, unsigned* ladder, 
                                     unsigned* wafer) {
    if (*tower > TKR_MAX_TOWER_ID) return false;
    if (tray == 0) return true;
    if (*tray > TKR_MAX_TRAY_ID)  return false;
    if (face == 0) return true;
    if (*face > TKR_MAX_FACE_ID)  return false;
    if (ladder == 0) return true;
    if (*ladder > TKR_MAX_LADDER_ID)  return false;
    if (wafer == 0) return true;
    return (*wafer <= TKR_MAX_WAFER_ID);
  }

  StatusCode 
  TkrInternalAlignCalib::getTrayAlign(unsigned towerId, unsigned trayId,
                                      Hep3Vector& disp, Hep3Vector& rot) {
    if (!checkLimits(&towerId, &trayId) ) return StatusCode::FAILURE;
    unsigned key = makeTrayKey(towerId, trayId);
    return getAlign(key, disp, rot);
  }

  StatusCode 
  TkrInternalAlignCalib::getFaceAlign(unsigned towerId, unsigned trayId, 
                                      unsigned faceId, Hep3Vector& disp, 
                                      Hep3Vector& rot) {

    if (!checkLimits(&towerId, &trayId, &faceId)) return StatusCode::FAILURE;
    unsigned key = makeFaceKey(towerId, trayId, faceId);
    return getAlign(key, disp, rot);
  }

  StatusCode 
  TkrInternalAlignCalib::getLadderAlign(unsigned towerId, unsigned trayId,
                                        unsigned faceId,unsigned ladderId,
                                        Hep3Vector& disp, Hep3Vector& rot)
  {
    if (!checkLimits(&towerId, &trayId, &faceId, &ladderId))
      return StatusCode::FAILURE;

    unsigned key = makeLadderKey(towerId, trayId, faceId, ladderId);
    return getAlign(key, disp, rot);
  }

  StatusCode 
  TkrInternalAlignCalib::getWaferAlign(unsigned towerId, unsigned trayId,
                                       unsigned faceId, unsigned ladderId,
                                       unsigned waferId, Hep3Vector& disp,
                                       Hep3Vector& rot) {
    if (!checkLimits(&towerId, &trayId, &faceId, &ladderId, &waferId))
      return StatusCode::FAILURE;

    unsigned key = makeWaferKey(towerId, trayId, faceId, ladderId, waferId);
    return getAlign(key, disp, rot);
  }

  /**
     Given key (which encodes item information, get associated constants)
     If no map entry for this key (or = null pointer) return 3-vectors
     of zeros.
  */
  StatusCode TkrInternalAlignCalib::getAlign(unsigned key, Hep3Vector& disp,
                                             Hep3Vector& rot)  {
    AlignMap::const_iterator pVal = m_items.find(key);
    if  (pVal != m_items.end()) {
      if (pVal->second != 0) {
        pVal->second->getAlign(disp, rot);
        return StatusCode::SUCCESS;
      }
    }
    assign3(Hep3Vector(0, 0, 0), disp);
    assign3(Hep3Vector(0, 0, 0), rot);

    return StatusCode::SUCCESS;
  }
  /**
     Given key (which encodes item information).  Make new item if need
     be; else replace existing.
  */
  StatusCode TkrInternalAlignCalib::putAlign(unsigned key, 
                                             const Hep3Vector& disp,
                                             const Hep3Vector& rot) {
    AlignMap::iterator pVal = m_items.find(key);
    bool found = false;
    if (pVal != m_items.end()) { //
      found = true;
      if (pVal->second !=0 ) {
        pVal->second->putAlign(disp, rot);
        return StatusCode::SUCCESS;
      }
    }
    TkrInternalAlign* pItem = new TkrInternalAlign(disp, rot);
    if (found) pVal->second = pItem;    // already was in map with null value
    else m_items[key] = pItem;
    return StatusCode::SUCCESS;
  }
  
  // Re-implemented from CalibBase
  StatusCode TkrInternalAlignCalib::update(CalibBase& other, MsgStream* log) {
    clean();
    TkrInternalAlignCalib& other1 = 
      dynamic_cast<TkrInternalAlignCalib&>(other);
    StatusCode sc = CalibBase::update(other1, log);
    if (sc != StatusCode::SUCCESS) return sc;

    m_items = other1.m_items;   // or do we need to do deep copy?

    return StatusCode::SUCCESS;               /* TO BE IMPLEMENTED */
  }




}
