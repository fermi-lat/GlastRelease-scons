#ifndef ACDTILEDIM_H
#define ACDTILEDIM_H

#include <vector>

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Geometry/Point3D.h"
#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"

class IGlastDetSvc;

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
  
  /// trivial destructor
  ~AcdTileDim() {;}
  
  /// direct access functions
  inline IGlastDetSvc& detSvc() const { return m_detSvc; }
  inline const idents::AcdId& acdId() const { return m_acdId; }
  inline const idents::VolumeIdentifier& volId() const { return m_volId; }
  inline const std::vector<double>& dim() const { return m_dim; }
  inline const HepPoint3D& tileCenter() const { return m_tileCenter; }
  inline const HepPoint3D* corner() const { return m_corners; }
  inline StatusCode statusCode() const { return m_sc; }

protected:
  
  /// this function access the detector service to get the geometry information
  StatusCode getVals();  

private:  

  /// The tile id
  const idents::AcdId             m_acdId;

  /// The volume id -> this is the key for the detector service
  const idents::VolumeIdentifier  m_volId;
  
  /// The detector service
  IGlastDetSvc&             m_detSvc;

  /// This show the status of the access to the detector service, should be checked before using data
  StatusCode                m_sc;

  /// the tile dimensions
  std::vector<double>       m_dim;

  /// the center of the tile (in global coordinates)
  HepPoint3D                m_tileCenter;

  /// the four corners of the tile (in global coordinates)
  HepPoint3D                m_corners[4];  

};
   

#endif
