#ifndef __AcdGeomMap_H
#define __AcdGeomMap_H 1

#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"

#include "AcdUtil/AcdRibbonDim.h"
#include "AcdUtil/AcdTileDim.h"
#include <map>
#include <set>

class IGlastDetSvc;

/**   
* @class AcdGeomMap
*
*
* $ $
*/


class AcdGeomMap {
	
public:
  
  AcdGeomMap(){;}
  
  virtual ~AcdGeomMap(){
    for ( std::map<idents::AcdId,AcdRibbonDim*>::iterator itrR = m_ribbonMap.begin();
	  itrR != m_ribbonMap.end(); itrR++ ) { delete itrR->second; }
    m_ribbonMap.clear();
    for ( std::map<idents::AcdId,AcdTileDim*>::iterator itrT = m_tileMap.begin();
	  itrT != m_tileMap.end(); itrT++ ) { delete itrT->second; }
    m_tileMap.clear();
  }
  
  /// @brief access
  const AcdRibbonDim* getRibbon(const idents::AcdId& id, const idents::VolumeIdentifier& volId, IGlastDetSvc &detSvc) {
    AcdRibbonDim* retVal(0);
    std::map<idents::AcdId,AcdRibbonDim*>::iterator itr = m_ribbonMap.find(id);
    if ( itr == m_ribbonMap.end() ) {
      retVal = new AcdRibbonDim(id,volId,detSvc);
      m_ribbonMap[id] = retVal;
      m_updated.insert(id);
    } else {
      retVal = itr->second;
      if ( m_updated.find(id) == m_updated.end() ) {
	retVal->update(detSvc);
      }
    }
    return retVal;
  }

  /// @brief access
  const AcdTileDim* getTile(const idents::AcdId& id, const idents::VolumeIdentifier& volId, IGlastDetSvc &detSvc) {
    AcdTileDim* retVal(0);
    std::map<idents::AcdId,AcdTileDim*>::iterator itr = m_tileMap.find(id);
    if ( itr == m_tileMap.end() ) {
      retVal = new AcdTileDim(id,volId,detSvc);
      m_tileMap[id] = retVal;
      m_updated.insert(id);
    } else {
      retVal = itr->second;
      if ( m_updated.find(id) == m_updated.end() ) {
	retVal->update(detSvc);
      }
    }
    return retVal;
  }

  void reset() {
    m_updated.clear();
  }
  
protected:
  
private:
  
  // stuff latched in from event
  std::map<idents::AcdId,AcdRibbonDim*> m_ribbonMap;
  std::map<idents::AcdId,AcdTileDim*>   m_tileMap;
  
  std::set<idents::AcdId> m_updated;

} ;

#endif
	
	
	
