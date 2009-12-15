#ifndef ACDRIBBONDIM_H
#define ACDRIBBONDIM_H

// stl headers
#include <vector>

// Gaudi, facilities, interfaces
#include "GaudiKernel/StatusCode.h"
#include "AcdUtil/IAcdGeometrySvc.h"

// Detector & geometry 
#include "idents/AcdId.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "geometry/Ray.h"


// A ribbon segment

class AcdRibbonSegment {
public:
  AcdRibbonSegment(float halfWidth, const HepPoint3D& start, const HepPoint3D& end, const HepTransform3D& trans);
  ~AcdRibbonSegment(){;}
  float          m_halfWidth;
  HepPoint3D     m_start;
  HepPoint3D     m_end;
  HepVector3D    m_vect;
  HepTransform3D m_trans;
  HepTransform3D m_invTrans;
};


/**
*  @class AcdRibbonDim
*
*  @brief  This class holds the dimension information about Acd Ribbons for use in the ACD reconstruction algorithms
*
*  In these algorithms the ribbons are as a series of rays.  Only the rays that actually cover gaps are used.  
*  This means that the small sections of ribbon that run perpendicular to the main ribbon direction in the
*  shingling region are ignored.
*
*  The ribbon width is taken as a constant everywhere.
*
*  In some cases we care about the total length along the ribbon.  That is measured from the center of the ribbon
*  values towards +X or +Y global are positive, towards -X or -Y global are negative
*  
*  \author Eric Charles
*
* $Header$
*/

class AcdRibbonDim {
  
public:
    
  typedef enum RibbonSide { MinusSide = -1,
			    Top = 0,
			    PlusSide = 1 };

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

  inline double halfLength() const { return m_halfLength; }
  inline StatusCode statusCode() const { return m_sc; }

  /**
   * @brief Convert a point for global to local coords
   *
   * @param global point in global frame
   * @param segment which ribbon segment
   * @param local point in local (ie segment) frame
   */
  void toLocal(const HepPoint3D& global, int segment, HepPoint3D& local) const;


  /**
   **/
  const AcdRibbonSegment* getSegment(unsigned iSeg) const {
    return m_segments[iSeg];
  }
  
  inline float getSegmentStart(unsigned iSeg) const {
    return m_segStarts[iSeg];
  }

  float getRibbonLength(unsigned iSeg, float segLength) const;


  inline unsigned int nSegments() const {
    return m_segments.size();
  }

  bool setEdgeRay(const int& side, HepPoint3D& start, HepVector3D& end) const;

  /**
   **/  
  void getSegmentsIndices(const HepVector3D& tkDir, bool upward, int& start, int& end, int& dir) const;

  /**
   * @brief Get the length along the ribbon given a volId and a local point
   *
   * @param volId which volume inside the ribbon
   * @param point the point in the volume frame
   * @param ribbonLength length along the ribbon (measured from the center)   
   * @param ribbonBin which bin to use for attenuation
   * @return true for success, false otherwise
   */
  bool getRibbonLengthAndBin(const idents::VolumeIdentifier &volId, const HepPoint3D& point,
			     double& ribbonLength, int& ribbonBin) const;

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

  /// size of the ribbons
  double                    m_halfLength;

  /// The index of the first top segment
  int                       m_topIdx;
  
  /// The index of the first +side segement
  int                       m_plusIdx;

  /// The various rays
  std::vector<AcdRibbonSegment*> m_segments;

  /// The starting points of the segments
  std::vector<float>        m_segStarts;
 

};
   

#endif
