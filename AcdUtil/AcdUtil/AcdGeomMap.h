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
*  @class AcdGeomMap
*
*  @brief  This class holds the geometrical information about the ACD elements.
*
*  The main function are getRibbon(...) and getTile(...) which return 
*  AcdRibbonDim and AcdTileDim object that encapsulate geometry information
*
* $Header$
*
*/


namespace AcdUtil {
  struct AcdVolumeAlignment {
  public:
    AcdVolumeAlignment(double x = 0., double y = 0., 
		       double sX = -1., double sY = -1.)
      :m_centerX(x),
       m_centerY(y),
       m_sizeX(sX),
       m_sizeY(sY){
    }    
    AcdVolumeAlignment(const AcdVolumeAlignment& other)
      :m_centerX(other.m_centerX),
       m_centerY(other.m_centerY),
       m_sizeX(other.m_sizeX),
       m_sizeY(other.m_sizeY){
    }
    AcdVolumeAlignment& operator=(const AcdVolumeAlignment& other){
      if ( this == &other ) return *this;
      m_centerX = other.centerX();
      m_centerY = other.centerY();
      m_sizeX = other.sizeX();
      m_sizeY = other.sizeY();
      return *this;
    }
    ~AcdVolumeAlignment(){;}    
    inline const double& centerX() const { return m_centerX; }
    inline const double& centerY() const { return m_centerY; }
    inline const double& sizeX() const { return m_sizeX; }
    inline const double& sizeY() const { return m_sizeY; }
  private:
    double m_centerX;
    double m_centerY;
    double m_sizeX;
    double m_sizeY;
  };
};


class AcdGeomMap {
	
public:
  
  /// Null c'tor
  AcdGeomMap()
    :m_acdGeomSvc(0){;}
  
  /// D'tor is just a simple cleanup
  virtual ~AcdGeomMap(){
    for ( std::map<idents::AcdId,AcdRibbonDim*>::iterator itrR = m_ribbonMap.begin();
	  itrR != m_ribbonMap.end(); itrR++ ) { delete itrR->second; }
    m_ribbonMap.clear();
    for ( std::map<idents::AcdId,AcdTileDim*>::iterator itrT = m_tileMap.begin();
	  itrT != m_tileMap.end(); itrT++ ) { delete itrT->second; }
    m_tileMap.clear();
  }
  
  /// @brief access
  const AcdRibbonDim* getRibbon(const idents::AcdId& id, IAcdGeometrySvc &detSvc) {
    AcdRibbonDim* retVal(0);    
    std::map<idents::AcdId,AcdRibbonDim*>::iterator itr = m_ribbonMap.find(id);
    if ( itr == m_ribbonMap.end() ) {
      //idents::AcdId& ncid = const_cast<idents::AcdId&>(id);
      retVal = new AcdRibbonDim(id,detSvc);
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
  const AcdTileDim* getTile(const idents::AcdId& id, IAcdGeometrySvc &detSvc) {
    AcdTileDim* retVal(0);
    std::map<idents::AcdId,AcdTileDim*>::iterator itr = m_tileMap.find(id);
    if ( itr == m_tileMap.end() ) {
      retVal = new AcdTileDim(id,detSvc);
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
  
  const AcdUtil::AcdVolumeAlignment& getAlignment(const idents::VolumeIdentifier& volId) const {
    static AcdUtil::AcdVolumeAlignment nullAlign;
    std::map<idents::VolumeIdentifier,AcdUtil::AcdVolumeAlignment>::const_iterator itr = m_alignMap.find(volId);
    return itr != m_alignMap.end() ? itr->second : nullAlign;
  }
  
  void putAlignment(const idents::VolumeIdentifier& volId, AcdUtil::AcdVolumeAlignment& acdAlign) {
    m_alignMap[volId] = acdAlign;
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

  std::map<idents::VolumeIdentifier,AcdUtil::AcdVolumeAlignment> m_alignMap;

  const IAcdGeometrySvc* m_acdGeomSvc;

} ;

#endif
	
	
	
