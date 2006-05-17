#ifndef ACDTILEDIM_H
#define ACDTILEDIM_H

#include <assert.h>

#include <vector>

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Geometry/Point3D.h"
#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CLHEP/Geometry/Transform3D.h"

/**
*  @class AcdTileDim
*
*  @brief This class holds the dimension information about Acd Tiles for use in the ACD reconstruction algorithms
*  In these algorithms the tiles are treated as thin, flat planes.  Therefore only the tile center and the four 
*  corners are used. 
*  
*  \author Eric Charles
*
* $Header$
*/

class AcdTileDim {

public:
    
  /// This function fills the corners using the information in the dimension vector and the center of the tile
  static StatusCode getCorners(const std::vector<double> &dim, const HepPoint3D &center, HepPoint3D *corner);

public:
    
  /// Constructor: takes the tile id, the volume id, and the detector service
  AcdTileDim(const idents::AcdId& acdId, const idents::VolumeIdentifier& volId, IGlastDetSvc &detSvc);

  /// Constructor: takes the tile id, both volume id, and the detector service
  AcdTileDim(const idents::AcdId& acdId, const idents::VolumeIdentifier& volIdMain, const idents::VolumeIdentifier& volIdOther, IGlastDetSvc &detSvc);
  
  /// trivial destructor
  ~AcdTileDim() {;}

  /// update (ie, when we get a new event)
  StatusCode update(IGlastDetSvc &detSvc) {
    m_detSvc = detSvc;
    m_sc = getVals();
    return m_sc;
  }

  /// direct access functions
  inline IGlastDetSvc& detSvc() const { return m_detSvc; }
  inline int nVol() const { return m_nVol; }
  inline const idents::AcdId& acdId() const { return m_acdId; }
  inline const idents::VolumeIdentifier& volId(int idx = 0) const { 
    assert(idx < m_nVol);
    return m_volId[idx]; 
  }
  inline const std::vector<double>& dim(int idx = 0) const { 
    assert(idx < m_nVol);
    return m_dim[idx]; 
  }
  inline const HepPoint3D& tileCenter(int idx =0) const { 
    assert(idx < m_nVol);
    return m_tileCenter[idx]; 
  }
  inline const HepPoint3D* corner(int idx = 0) const { 
    assert(idx < m_nVol);
    return m_corners[idx]; 
  }
  inline StatusCode statusCode() const { return m_sc; }

  inline int sharedEdge(int idx) const { 
    assert(idx < m_nVol);
    return m_shared[idx]; 
  }

  inline float sharedWidth(int idx) const { 
    assert(idx < m_nVol);
    return m_sharedWidth[idx]; 
  }

  inline int face(int idx) const {
    assert(idx < m_nVol);
    return (m_volId[idx])[1];    
  }

  void toLocal(const HepPoint3D& global, HepPoint3D& local, int idx = 0);

protected:
  
  /// this function access the detector service to get the geometry information
  StatusCode getVals();  

private:  

  /// The tile id
  const idents::AcdId             m_acdId;

  /// The number of volumes in this tile
  const int                 m_nVol;

  /// The volume id -> this is the key for the detector service
  idents::VolumeIdentifier  m_volId[2];
  
  /// The detector service
  IGlastDetSvc&             m_detSvc;

  /// This show the status of the access to the detector service, should be checked before using data
  StatusCode                m_sc;

  /// the tile dimensions
  std::vector<double>       m_dim[2];

  /// the center of the tile (in global coordinates)
  HepPoint3D                m_tileCenter[2];

  /// the four corners of the tile (in global coordinates)
  HepPoint3D                m_corners[2][4];  

  /// the transformations to local coords
  HepTransform3D            m_transform[2];

  /// which (if any) edges are shared between tiles:
  int                       m_shared[2];
  
  /// width extra of extra volume in curved tiles
  float                     m_sharedWidth[2];

};
   

#endif
