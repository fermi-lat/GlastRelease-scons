#ifndef ACDRECONFUNCS_H
#define ACDRECONFUNCS_H

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Geometry/Transform3D.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "geometry/Ray.h"

#include "./AcdReconStruct.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"

class AcdTileDim;
class AcdTkrParams;
class AcdRibbonDim;
class IPropagator;

namespace AcdUtil {
  class RayDoca;
}

namespace AcdRecon {

  // A couple of constants
  const int ribbonX = 5;
  const double maxDocaValue = 2000.;

  /**
   * @brief POCA Between a track and a point
   *
   * @param track the track projection data
   * @param point the point in question
   * @param arcLength at which the POCA occurs
   * @param doca distance of closest approach
   * @param poca point of closest approach
   */
  void pointPoca(const AcdRecon::TrackData& track, const HepPoint3D& point,
		 double& arcLength, double& doca, HepPoint3D& poca);

  /**
   * @brief POINT where a track crosses a plane of constant X,Y,Z)

   *
   * @param track the track projection data
   * @param plane point at the center of the plane, in GLOBAL coordinates
   * @param face orientation of the place
   * @param arcLength at which the intersection occurs
   * @param hitPoint of intersection
   */
  void crossesPlane(const AcdRecon::TrackData& track, const HepPoint3D& plane, int face, 
		    double& arcLength, HepPoint3D& hitPoint);  
  
  /**
   * @brief POINT where a track crosses a plane defined by a HEP transformation
   *
   * @param track the track projection data
   * @param plane transformation to the center of the plane, include rotation
   * @param arcLength at which the intersection occurs
   * @param hitPoint of intersection
   */
  void crossesPlane(const AcdRecon::TrackData& track, const HepTransform3D& plane, 
		    double& arcLength, HepPoint3D& hitPoint);  

  /**
   * @brief POCA Between a track and a ray defined by two points
   *
   * @param track the track projection data
   * @param c1 first corner
   * @param c2 second corner
   * @param rayDoca data about the POCA
   * @param edgeLen length along the ray between the points to the POCA
   */
  void rayDoca(const AcdRecon::TrackData& track, const HepPoint3D& c1, const HepPoint3D& c2,
	       AcdUtil::RayDoca& rayDoca, double& edgeLen);
  
  /**
   * @brief POCA Between a track and a ray, 
   *        returns DOCA to point at end of ray if DOCA occurs outside ray edges
   *
   * @param track the track projection data
   * @param ray the ray in question
   * @param arcLength at which the POCA occurs
   * @param rayLength length along ray at which the POCA occurs
   * @param dist distance of closest approach
   * @param x POCA of track to ray.  (ie, a closest point on track to ray)
   * @param v vector from POCA (x-v must be on ray)
   * @param region 
   */
  // 
  // 
  void rayDoca_withCorner(const AcdRecon::TrackData& track, const Ray& ray,
			  double& arcLength, double& rayLength, double& dist, Point& x, Vector& v, int& region);

  /**
   * @brief POINT where a track crosses the plane of a tile
   *
   * @param track the track projection data
   * @param tile AcdTileDim object with tile geometry
   * @param data AcdRecon::PocaData object that encapsulates intersection data
   *   This function fills:
   * 
   */
  void tilePlane(const AcdRecon::TrackData& track, const AcdTileDim& tile,
		 AcdRecon::PocaData& data);

  /**
   * @brief Transform tile plane crossing point into active distance
   *
   * @param tile AcdTileDim object with tile geometry
   * @param iVol normally 0, 1 for bent pieces
   * @param globalPoint point in global frame
   * @param localPoint point in local frame
   * @param activeX distance to closest edge in local X, > 0 -> inside, < 0 -> outside
   * @param activeY distance to closest edge in local Y, > 0 -> inside, < 0 -> outside
   */
  void tilePlaneActiveDistance(const AcdTileDim& tile, int iVol, const HepPoint3D& globalPoint,
			       HepPoint3D& localPoint, double& activeX, double& activeY);

  /**
   * @brief POCA Between a track and edge of a tile (when the track goes inside the tile)
   *
   * @param track the track projection data
   * @param tile AcdTileDim object with tile geometry
   * @param arcLength at which the POCA occurs
   * @param dist distance of closest approach
   * @param x POCA of track to ray.  (ie, a closest point on track to ray)
   * @param v vector from POCA (x-v must be on ray)
   * @param region 
   */
  void tileEdgePoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
		    double& arcLength, double& dist, Point& x, Vector& v, int& region);

  /**
   * @brief POCA Between a track and edge of a tile (when the track goes outside the tile)
   *
   * @param track the track projection data
   * @param tile AcdTileDim object with tile geometry
   * @param arcLength at which the POCA occurs
   * @param dist distance of closest approach
   * @param x POCA of track to ray.  (ie, a closest point on track to ray)
   * @param v vector from POCA (x-v must be on ray)
   * @param region 
   */
  void tileEdgeCornerPoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
			  double& arcLength, double& dist, Point& x, Vector& v, int& region);  

  /**
   * @brief 
   *
   * @param track the track projection data
   * @param ribbon AcdRibbonDim object with ribbon geometry
   * @param data AcdRecon::PocaData object that encapsulates intersection data
   *   This function fills:
   */
  // POINT where a track crosses the plane of a ribbon
  void ribbonPlane(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon,
		   AcdRecon::PocaData& data);


  /**
   * @brief  POCA between a track and a ribbon 
   * (includes all the ribbon segments in the ribbon direction, but not the small perpindicular segements
   *
   * @param track the track projection data
   * @param ribbon AcdRibbonDim object with ribbon geometry
   * @param arcLength at which the POCA occurs
   * @param rayLength length along ribbon at which the POCA occurs
   * @param dist distance of closest approach
   * @param x POCA of track to ray.  (ie, a closest point on track to ray)
   * @param v vector from POCA (x-v must be on ray)
   * @param region 
   */
  void ribbonPoca(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon,
		  double& arcLength, double& rayLength, double& dist, Point& x, Vector& v, int& region);  

  /**
   * @brief initiliaze the kalman propagator
   *
   * @param prop the propagator
   * @param aTrack the full TKR track representation
   * @param trackData the minimal ACD track representation
   * @param maxArcLength furthers arcLength we want to extend the track to
   */
  void startPropagator(IPropagator& prop, const Event::TkrTrack& aTrack, const AcdRecon::TrackData& trackData,
		       const double& maxArcLength);
  
  /**
   * @brief 
   *
   * @param prop the propagator
   * @param trackDatat the minimal ACD track representation
   * @param arcLength current arclength of the track
   * @param next_params full TKR rep of parameters at next step in Kalman fit
   * @param paramsAtArcLength minimal ACD rep of parameters at next step
   */
  // run the propagtor out to a specified arclength
  void propagateToArcLength(IPropagator& prop,
			    const AcdRecon::TrackData& trackData, const double& arcLength,
			    Event::TkrTrackParams& next_params,
			    AcdTkrParams& paramsAtArcLength);

  /**
   * @brief Project error ellipsoid onto a plane
   *
   * @param paramsAtArcLength minimal ACD rep of track parameters parameters 
   * @param localFrameVectors vectors pointing along local X and Y axis, in global frame
   * @param covAtPlane covarience matrix at that place
   */
  void projectErrorToPlane(const AcdTkrParams& paramsAtArcLength, const CLHEP::HepMatrix& localFrameVectors,
			   CLHEP::HepSymMatrix& covAtPlane);

  /**
   * @brief  Project error ellipsoid onto a line
   *
   * @param paramsAtArcLength minimal ACD rep of track parameters parameters 
   * @param pocaVector vector from POCA to closest point on tile
   * @param pocaError error projected along pocaVector
   */
  void projectErrorToPocaVector(const AcdTkrParams& paramsAtArcLength, const Vector& pocaVector, 
				double& pocaError);

  /**
   * @brief Project where a track exits the ACD volume (aka the rectangular solid defined by the ACD)
   * Note that this does either the upgoing or downgoing intersection, depending on which ray
   * trackData is defined as
   *
   * @param trackData the track projection data
   * @param acdVol simplified ACD geometry (just a box)
   * @param data AcdRecon::ExitData object that encapsulates intersection data
   *   This function fills:
   */
  bool exitsLat(const AcdRecon::TrackData& trackData,
		const AcdRecon::AcdVolume& acdVol,
		AcdRecon::ExitData& data);
  
  /**
   * @brief Project where a track enters the ACD volume (aka the rectangular solid defined by the ACD)
   * This is mainly used from MC particles, since they originate outside the ACD
   * This returns only the first point where track pierces ACD
   *
   * @param trackData the track projection data
   * @param acdVol simplified ACD geometry (just a box)
   * @param data AcdRecon::ExitData object that encapsulates intersection data
   *   This function fills:
   */
  bool entersLat(const AcdRecon::TrackData& trackData,
		 const AcdRecon::AcdVolume& acdVol,
		 AcdRecon::ExitData& data); 
  
  /**
   * @brief Project where a track enters the CAL volume (or really, where it crosses a plane of constant
   * Z at the to of the CAL)
   *
   * @param trackData the track projection data
   * @param calZDist z value where the CAL starts
   * @param entryPoint point where track enters CAL
   * @param entryVector vector along track where track enters CAL
   * @param region
   */
  void entersCal(const AcdRecon::TrackData& trackData, const double& calZDist, 
		 Point& entryPoint, Vector& entryVector, int& region);
 
  /**
   * @brief 
   *
   * @param tile AcdTileDim object with tile geometry
   * @param entryPoint point where track enters CAL
   * @param entryVector vector along track where track enters CAL
   * @param solidAngle solid angle subtended by tile, w.r.t. entryPoint
   * @param meanAngle mean angle between entryVector and entryPoint w.r.t. tile surface
   * @param pathLength 
   */
  void splashVariables(const AcdTileDim& tile, const Point& entryPoint, const Vector& entryVector,
		       double& solidAngle, double& meanAngle, double& pathLength);
 
}

#endif
