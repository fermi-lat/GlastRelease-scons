#ifndef ACDTILEDIM_H
#define ACDTILEDIM_H

// stl headers
#include <vector>
#include <assert.h>

// Gaudi, facilities, interfaces
#include "GaudiKernel/StatusCode.h"
#include "AcdUtil/IAcdGeometrySvc.h"

// Detector & geometry 
#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Matrix/Matrix.h"


/**
*  @class AcdTileSection
*
*  @brief This class holds the dimension information a part of an ACD Tile
*  
*  In these algorithms the tiles are treated as 1 or 2 thin, flat planes.  AcdTileSection describes one of those planes.
*
*  This class has 2 represetations of the tile:
*
*    1) 4 corners in global coodinates
*    2) Center in global cordinates, rotation matrix, and dimensions in local coordinates.
*
*
*  Since some of the tiles slightly trapezoidal, with the width as we move along the local Y frame
*  we provide functions to estimate the active distance and the width at any point near the tile.
*
*  
*  
*  \author Eric Charles
*
* $Header$
**/


class AcdTileSection {
public:
  AcdTileSection();
  ~AcdTileSection(){;}

  double widthAtLocalY(double localY) const;
    
  void activeDistance(const HepPoint3D& localPoint, int iVol, double& activeX, double& activeY) const;

  std::vector<double> m_dim;  
  bool                m_trapezoid;
  float               m_xmean;
  float               m_xdiff;
  int                 m_shared;
  float               m_sharedWidth;  
  HepPoint3D          m_center;
  HepPoint3D          m_corners[4];
  HepTransform3D      m_trans;
  HepTransform3D      m_invTrans; 
 
};


/**
*  @class AcdTileDim
*
*  @brief This class holds the dimension information about Acd Tiles for use in the ACD reconstruction algorithms
*  
*  In these algorithms the tiles are treated as thin, flat planes.  Therefore only the tile center and the four 
*  corners are used. 
*
*  Most of the tiles are a single piece, however 10 tiles, (0-4 and 40-44) are made from two pieces, one on the top
*  of the ACD and one that runs down the side.  This makes things quite a bit more complicated.
*  
*  \author Eric Charles
*
* $Header$
**/

class AcdTileDim {

public:
    
  /// Constructor: takes the tile id, the volume id, and the detector service
  AcdTileDim(const idents::AcdId& acdId, IAcdGeometrySvc& acdGeomSvc);
  
  /// trivial destructor
  virtual ~AcdTileDim();

  /// update (ie, when we get a new event)
  StatusCode update(IAcdGeometrySvc& acdGeomSvc) {
    m_acdGeomSvc = acdGeomSvc;
    m_sc = getVals();
    return m_sc;
  }

  StatusCode latchGaps();

  /// direct access functions
  inline IAcdGeometrySvc& acdGeomSvc() const { return m_acdGeomSvc; }
  inline const idents::AcdId& acdId() const { return m_acdId; }
  inline StatusCode statusCode() const { return m_sc; }

  inline int nVol() const { return m_sections.size(); }

  inline const AcdTileSection* getSection(int iVol) const {
    return m_sections[iVol];
  }

  inline const std::vector< HepPoint3D > screwHoles() const {
    return m_screwHoles;
  }

  float gapAtLocalY(int side, double localY) const;
    
  float gapSize(const HepPoint3D& localPoint, int iVol) const;

  inline const HepPoint3D& tileCenter(unsigned int iVol) const {
    assert(iVol < m_sections.size());
    return m_sections[iVol]->m_center;
  }

  inline const HepTransform3D& transform(unsigned int iVol) const {
    assert(iVol < m_sections.size());
    return m_sections[iVol]->m_trans;
  }

  inline const HepPoint3D* corner(unsigned int iVol) const {
    assert(iVol < m_sections.size());
    return m_sections[iVol]->m_corners;
  }

  inline void toLocal(const HepPoint3D& global, HepPoint3D& localPoint, unsigned int iVol) const {
    assert(iVol < m_sections.size());
    localPoint = m_sections[iVol]->m_trans * global;
  }

  
protected:
  
  /// this function access the detector service to get the geometry information
  StatusCode getVals();  

private:  

  /// The tile id
  const idents::AcdId       m_acdId;

  /// ACD Geom Service
  IAcdGeometrySvc&          m_acdGeomSvc;

  /// This show the status of the access to the detector service, should be checked before using data
  StatusCode                m_sc;

  /// The Sections
  std::vector<AcdTileSection*> m_sections;

  /// Positions of the screw holes, in local coords
  std::vector< HepPoint3D > m_screwHoles;

  bool                m_doneGaps;

  float               m_minusXGap[2];
  float               m_plusXGap[2];
 
};
   

#endif
