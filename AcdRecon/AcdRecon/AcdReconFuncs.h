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

  //
  void pointPoca(const AcdRecon::TrackData& track, const HepPoint3D& point,
		 double& arcLength, double& doca, HepPoint3D& poca);

  void crossesPlane(const AcdRecon::TrackData& track, const HepPoint3D& plane, int face, 
		    double& arcLength, double& localX, double& localY, HepPoint3D& hitPoint);  

  void rayDoca(const AcdRecon::TrackData& track, const HepPoint3D& c1, const HepPoint3D& c2,
	       RayDoca& rayDoca, double& edgeLen);

  void tilePlane(const AcdRecon::TrackData& track, const AcdTileDim& tile,
		 double& arcLength, double& localX, double& localY, 
		 double& activeX, double& activeY, double& active2D, HepPoint3D& hitPoint);  

  void tileEdgePoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
		    double& arcLength, double& dist, Point& x, Vector& v, int& region);

  void tileEdgeCornerPoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
			  double& arcLength, double& dist, Point& x, Vector& v, int& region);  

  void ribbonPlane(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon, 
		  double& arcLength, double& dist, HepPoint3D& x);

  void ribbonPoca(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon,
		  double& arcLength, double& dist, Point& x, Vector& v, int& region);

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
