#ifndef ACDRECONFUNCS_H
#define ACDRECONFUNCS_H

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Geometry/Transform3D.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "./AcdReconStruct.h"
#include "./RayDoca.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"

class AcdTileDim;
class AcdTkrParams;
class AcdRibbonDim;
class IPropagator;

namespace AcdRecon {

  // A couple of constants
  const int ribbonX = 5;
  const double maxDocaValue = 2000.;

  // POCA Between a track and a point
  void pointPoca(const AcdRecon::TrackData& track, const HepPoint3D& point,
		 double& arcLength, double& doca, HepPoint3D& poca);

  // POINT where a track crosses a plane of constant X,Y,Z)
  void crossesPlane(const AcdRecon::TrackData& track, const HepPoint3D& plane, int face, 
		    double& arcLength, HepPoint3D& hitPoint);  
  
  // POINT where a track crosses a plane defined by a HEP transformation
  void crossesPlane(const AcdRecon::TrackData& track, const HepTransform3D& plane, 
		    double& arcLength, HepPoint3D& hitPoint);  

  // POCA Between a track and a ray defined by two points
  // edgeLen is the length along the ray between the points
  void rayDoca(const AcdRecon::TrackData& track, const HepPoint3D& c1, const HepPoint3D& c2,
	       RayDoca& rayDoca, double& edgeLen);
  
  // POCA Between a track and a ray
  // returns DOCA to point at end of ray if DOCA occurs outside ray edges
  void rayDoca_withCorner(const AcdRecon::TrackData& track, const Ray& ray,
			  double& arcLength, double& rayLength, double& dist, Point& x, Vector& v, int& region);

  // POINT where a track crosses the plane of a tile
  void tilePlane(const AcdRecon::TrackData& track, const AcdTileDim& tile,
		 AcdRecon::PocaData& data);

  // Transform tile plane crossing point into active distance
  void tilePlaneActiveDistance(const AcdTileDim& tile, int iVol, const HepPoint3D& globalPoint,
			       HepPoint3D& localPoint, double& activeX, double& activeY);

  // POCA Between a track and edge of a tile (when the track goes inside the tile)
  void tileEdgePoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
		    double& arcLength, double& dist, Point& x, Vector& v, int& region);

  // POCA Between a track and the edges and corners of a tile (when the track goes outside the tile) 
  void tileEdgeCornerPoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
			  double& arcLength, double& dist, Point& x, Vector& v, int& region);  

  // POINT where a track crosses the plane of a ribbon
  void ribbonPlane(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon,
		   AcdRecon::PocaData& data);


  // POCA between a track and a ribbon 
  // (includes all the ribbon segments in the ribbon direction, but not the small perpindicular segements
  void ribbonPoca(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon,
		  double& arcLength, double& rayLength, double& dist, Point& x, Vector& v, int& region);  

  // initiliaze the kalman propagator
  void startPropagator(IPropagator& prop, const Event::TkrTrack& aTrack, const AcdRecon::TrackData& trackData,
		       const double& maxArcLength);
  
  // run the propagtor out to a specified arclength
  void propagateToArcLength(IPropagator& prop,
			    const AcdRecon::TrackData& trackData, const double& arcLength,
			    Event::TkrTrackParams& next_params,
			    AcdTkrParams& paramsAtArcLength);

  // Project error ellipsoid onto a plane
  void projectErrorToPlane(const AcdTkrParams& paramsAtArcLength, const CLHEP::HepMatrix& localFrameVectors,
			   CLHEP::HepSymMatrix& covAtPlane);

  // Project error ellipsoid onto a line
  void projectErrorToPocaVector(const AcdTkrParams& paramsAtArcLength, const Vector& pocaVector, 
				double& pocaError);

  // Project where a track exits the ACD volume (aka the rectangular solid defined by the ACD)
  // Note that this does either the upgoing or downgoing intersection, depending on which ray
  // trackData is defined as
  bool exitsLat(const AcdRecon::TrackData& trackData,
		const AcdRecon::AcdVolume& acdVol,
		AcdRecon::ExitData& data);
  
  // Project where a track enters the ACD volume (aka the rectangular solid defined by the ACD)
  // This is mainly used from MC particles, since they originate outside the ACD
  // This returns only the first point where track pierces ACD
  bool entersLat(const AcdRecon::TrackData& trackData,
		 const AcdRecon::AcdVolume& acdVol,
		 AcdRecon::ExitData& data); 
  
  void entersCal(const AcdRecon::TrackData& trackData, const double& calZDist, 
		 Point& entryPoint, Vector& entryVector, int& region);
 
  void splashVariables(const AcdTileDim& tile, const Point& entryPoint, const Vector& entryVector,
		       double& solidAngle, double& meanAngle, double& pathLength);
 
}

#endif
