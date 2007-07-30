#ifndef ACDRECONFUNCS_H
#define ACDRECONFUNCS_H

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Matrix/Matrix.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdRibbonDim.h"

#include "./AcdReconStruct.h"
#include "./RayDoca.h"

#include "Event/Recon/TkrRecon/TkrTrackParams.h"


namespace AcdRecon {

  // A couple of constats
  const int ribbonX = 5;
  const double maxDocaValue = 2000.;

  // POCA Between a track and a point
  void pointPoca(const AcdRecon::TrackData& track, const HepPoint3D& point,
		 double& arcLength, double& doca, HepPoint3D& poca);

  // POINT where a track crosses a plane of constant X,Y,Z)
  void crossesPlane(const AcdRecon::TrackData& track, const HepPoint3D& plane, int face, 
		    double& arcLength, double& localX, double& localY, HepPoint3D& hitPoint);  

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
		 double& arcLength, double& localX, double& localY, 
		 double& activeX, double& activeY, double& active2D, HepPoint3D& hitPoint);  

  // POCA Between a track and edge of a tile (when the track goes inside the tile)
  void tileEdgePoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
		    double& arcLength, double& dist, Point& x, Vector& v, int& region);

  // POCA Between a track and the edges and corners of a tile (when the track goes outside the tile) 
  void tileEdgeCornerPoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
			  double& arcLength, double& dist, Point& x, Vector& v, int& region);  

  // POINT where a track crosses the plane of a ribbon
  void ribbonPlane(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon, 
		   double& arcLength, double& dist, HepPoint3D& x);

  // POCA between a track and a ribbon 
  // (includes all the ribbon segments in the ribbon direction, but not the small perpindicular segements
  void ribbonPoca(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon,
		  double& arcLength, double& rayLength, double& dist, Point& x, Vector& v, int& region);  

  // Point where a track g
  void projectToPlane(const AcdRecon::TrackData& trackData, const Event::TkrTrackParams& trackParam, 
		      int face, const HepPoint3D& plane,
		      AcdRecon::PocaData& data);

  void errorAtXPlane(double delta, const Event::TkrTrackParams& track, CLHEP::HepMatrix& covAtPlane);
  
  void errorAtYPlane(double delta, const Event::TkrTrackParams& track, CLHEP::HepMatrix& covAtPlane);
  
  void errorAtZPlane(double delta, const Event::TkrTrackParams& track, CLHEP::HepMatrix& covAtPlane);

  void projectErrorAtPoca(const AcdRecon::TrackData& trackData, const Event::TkrTrackParams& trackParams,
			  const Point& poca, const Vector& pocaVector, double& pocaError);

  void entersCal(const AcdRecon::TrackData& trackData, const double& calZDist, 
		 Point& entryPoint, Vector& entryVector, int& region);
 
  void splashVariables(const AcdTileDim& tile, const Point& entryPoint, const Vector& entryVector,
		       double& solidAngle, double& meanAngle, double& pathLength);
 
}

#endif
