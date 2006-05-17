#ifndef __AcdGeomMap_H
#define __AcdGeomMap_H 1

#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"

#include "AcdUtil/AcdRibbonDim.h"
#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/IAcdGeometrySvc.h"
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
  
  AcdGeomMap()
    :m_acdGeomSvc(0){;}
  
  virtual ~AcdGeomMap(){
    for ( std::map<idents::AcdId,AcdRibbonDim*>::iterator itrR = m_ribbonMap.begin();
	  itrR != m_ribbonMap.end(); itrR++ ) { delete itrR->second; }
    m_ribbonMap.clear();
    for ( std::map<idents::AcdId,AcdTileDim*>::iterator itrT = m_tileMap.begin();
	  itrT != m_tileMap.end(); itrT++ ) { delete itrT->second; }
    m_tileMap.clear();
  }
  
  /// @brief access
  const AcdRibbonDim* getRibbon(const idents::AcdId& id, IGlastDetSvc &detSvc) {
    AcdRibbonDim* retVal(0);    
    std::map<idents::AcdId,AcdRibbonDim*>::iterator itr = m_ribbonMap.find(id);
    if ( itr == m_ribbonMap.end() ) {
      const idents::VolumeIdentifier volId = id.volId();
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
  const AcdTileDim* getTile(const idents::AcdId& id, IGlastDetSvc &detSvc) {
    AcdTileDim* retVal(0);
    std::map<idents::AcdId,AcdTileDim*>::iterator itr = m_tileMap.find(id);
    if ( itr == m_tileMap.end() ) {
      const std::map<idents::AcdId, int>::const_iterator itrCount = m_acdGeomSvc->getAcdIdVolCountCol().find(id);
      int nVol = itrCount->second;
      idents::VolumeIdentifier volId = id.volId();
      idents::VolumeIdentifier volIdSide;
      switch (nVol){
      case 1:	
	retVal = new AcdTileDim(id,volId,detSvc);
	break;
      case 2:
	volIdSide = id.volId(true);
	retVal = new AcdTileDim(id,volId,volIdSide,detSvc);
	break;
      }
      if ( retVal == 0 ) return retVal;
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
  
  void setAcdGeomSvc(const IAcdGeometrySvc& svc){
    m_acdGeomSvc = &svc;
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

  const IAcdGeometrySvc* m_acdGeomSvc;

} ;

#endif
	
	
	
