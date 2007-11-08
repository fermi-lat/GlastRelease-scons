#ifndef ACDRIBBONDIM_H
#define ACDRIBBONDIM_H

#include <vector>

#include "GaudiKernel/StatusCode.h"
#include "idents/AcdId.h"
#include "CLHEP/Geometry/Transform3D.h"

#include "AcdUtil/IAcdGeometrySvc.h"
#include "geometry/Ray.h"


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
    
  /// Constructor: takes the tile id, the volume id, and the detector service
  AcdRibbonDim(const idents::AcdId& acdId, IAcdGeometrySvc& acdGeomSvc);
  
  /// trivial destructor  
  ~AcdRibbonDim() {;}

  /// update (ie, when we get a new event)
  StatusCode update(IAcdGeometrySvc& acdGeomSvc) {
    bool newVals(true);
    m_acdGeomSvc = acdGeomSvc;
    if ( newVals ) {
      m_sc = getVals();
    }
    return m_sc;
  }

  /// direct access functions
  inline IAcdGeometrySvc& acdGeomSvc() const { return m_acdGeomSvc; }
  inline const idents::AcdId& acdId() const { return m_acdId; }

  inline double halfWidth() const { return m_halfWidth; }
  inline double halfLength() const { return m_halfLength; }
  inline StatusCode statusCode() const { return m_sc; }
  inline const std::vector<Ray>& minusSideRays() const { return m_minusSideRays; }     
  inline const std::vector<Ray>& topRays() const { return m_topRays; }      
  inline const std::vector<Ray>& plusSideRays() const { return m_plusSideRays; }        

  bool setEdgeRay(int iSeg, HepPoint3D& start, HepVector3D& vector) const;

  void toLocal(const HepPoint3D& global, int segment, HepPoint3D& local) const;

  const HepTransform3D& transform(int iVol) const;

protected:

  /// this function access the detector service to get the geometry information
  StatusCode getVals();  

  /// this function adds up the length of the rays the old-fashioned way, cause that is safer
  double calcRibbonLength();

private:  

  /// The tile id
  const idents::AcdId       m_acdId;

  /// ACD Geom Service      
  IAcdGeometrySvc&          m_acdGeomSvc;

  /// This show the status of the access to the detector service, should be checked before using data
  StatusCode                m_sc;

  /// the transformations to local coords
  HepTransform3D            m_minusSideTransform;  
  HepTransform3D            m_topTransform;
  HepTransform3D            m_plusSideTransform;

  /// size of the ribbons
  double                    m_halfWidth;  
  double                    m_halfLength;

  /// The various rays
  std::vector<Ray>          m_minusSideRays;
  std::vector<Ray>          m_topRays;
  std::vector<Ray>          m_plusSideRays;         
 
};
   

#endif
