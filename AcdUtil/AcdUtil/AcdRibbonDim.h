#ifndef ACDRIBBONDIM_H
#define ACDRIBBONDIM_H

#include <vector>

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Geometry/Point3D.h"
#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

class IGlastDetSvc;

#ifndef HepPoint3D
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

/**
*  @class AcdRibbonDim
*
*
*  @brief  This class holds the dimension information about Acd Ribbons for use in the ACD reconstruction algorithms
*  In these algorithms the ribbons are treated 3 lines.  Therefore only the line start and end point are used. 
*  
*  \author Eric Charles
*
* $Header$
*/

class AcdRibbonDim {

public:
    
  /// This function gets information about the ribbons from the detector service
  StatusCode getDetectorDimensions( int isegment, 
				    const idents::VolumeIdentifier& volId, IGlastDetSvc &detSvc,
				    std::vector<double>& dim, HepPoint3D& center);
  
public:
    
  /// Constructor: takes the tile id, the volume id, and the detector service
  AcdRibbonDim(const idents::AcdId& acdId, const idents::VolumeIdentifier& volId, IGlastDetSvc &detSvc);
  
  /// trivial destructor  
  ~AcdRibbonDim() {;}

  /// update (ie, when we get a new event)
  StatusCode update(IGlastDetSvc &detSvc) {
    m_detSvc = detSvc;
    m_sc = getVals();
    return m_sc;
  }

  /// direct access functions
  inline IGlastDetSvc& detSvc() const { return m_detSvc; }
  inline const idents::AcdId& acdId() const { return m_acdId; }
  inline const idents::VolumeIdentifier& volId() const { return m_volId; }
  inline const idents::VolumeIdentifier* segId() const { return m_segId; }
  inline const HepPoint3D* ribbonStart() const { return m_start; }
  inline const HepPoint3D* ribbonEnd() const { return m_end; }
  inline const double* halfWidth() const { return m_halfWidth; }
  inline StatusCode statusCode() const { return m_sc; }

  void toLocal(const HepPoint3D& global, int segment, HepPoint3D& local);

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

  /// The volume ids for the segments -> these are more keys for the detector service
  idents::VolumeIdentifier  m_segId[3];

  /// the start points of the segments (in global coordinates)
  HepPoint3D                m_start[3];
  
  /// the end points of the segments (in global coordinates)
  HepPoint3D                m_end[3];

  /// the transformations to local coords
  HepTransform3D            m_transform[3];

  /// size of the ribbons
  double                    m_halfWidth[3];

};
   

#endif
