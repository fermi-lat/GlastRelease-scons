#ifndef ACDRECONFUNCSV2_H
#define ACDRECONFUNCSV2_H

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
#include "Event/Recon/TkrRecon/TkrTrackParams.h"
#include "Event/Recon/CalRecon/CalParams.h"


class AcdTileDim;
class AcdTkrParams;
class AcdRibbonDim;
class IPropagator;

namespace AcdUtil {
  class RayDoca;
}

namespace AcdRecon {

  namespace ReconFunctions {

    using namespace CLHEP;
    typedef HepGeom::Point3D<double> HepPoint3D;

    // A couple of constants
    const double maxDocaValue = 2000.;
    
    /**
     * @brief fill the 5 x 4 derivative matrix to get from Tkr covariance form to Acd covariance form
     *
     * @param trackDir the directional cosines
     * @param cov the new covariance matrix
     **/
    void fillTkrToAcdCovTranslation(const HepVector3D& dir, CLHEP::HepMatrix& covTrans);

    /**
     * @brief fill the 4 x 5 derivative matrix to get from Acd covariance form to Tkr covariance form
     *
     * @param dir the directional cosines
     * @param covTrans the new covariance matrix
     **/
    void fillAcdToTkrCovTranslation(const HepVector3D& dir, CLHEP::HepMatrix& covTrans);

    /**
     * @brief fill the 4 x 4 covariance matrix from TrkTrackParams
     *
     * @param tkrParams the directional cosines
     * @param cov the  matrix
     **/
    void fillCovMatrixFromTkr(const Event::TkrTrackParams& trackParams, CLHEP::HepSymMatrix& cov);
    
    /**
     * @brief Convert from TrkTrackParams to AcdRecon::TrackData
     *
     * @param trackParams the tracker params
     * @param zRef value of z at ref point
     * @param energy at starting point
     * @param index of track
     * @param up true if projection is upgoing
     * @param acdParams the same, in AcdRecon::TrackData rep
     **/
    void convertToAcdRep(const Event::TkrTrackParams& trackParams,
			 double zRef,
			 AcdRecon::TrackData& acdParams);

    /**
     * @brief Convert from AcdRecon::TrackData to TrkTrackParams
     *
     * @param acdParams the AcdRecon::TrackData rep
     * @param trackParams the tracker params
     **/
    void convertToTkrRep(const AcdRecon::TrackData& acdParams,
			 Event::TkrTrackParams& trackParams);


    /**
     * @brief Convert from CalParams to AcdRecon::TrackData
     *
     * @param calParams the cal cluster params
     * @param zRef value of z at ref point
     * @param acdParams the same, in AcdRecon::TrackData rep
     **/
    void convertToAcdRep(Event::CalParams& calParams,
			 AcdRecon::TrackData& acdParams);
    
    /**
     * @brief POCA Between a track and a point
     *
     * @param track the track projection data
     * @param point the point in question
     * @param arcLength at which the POCA occurs
     * @param poca point of closest approach
     * @param voca vector of closest approach
     * @return true for success, false if poca behind head of track
     **/
    bool pointPoca(const AcdRecon::TrackData& track, const HepPoint3D& point,
		   double& arcLength, HepPoint3D& poca, HepVector3D& voca );
    
    
    /**
     * @brief Error on the POCA Between a track and a point
     *
     * @param track the track projection data
     * @param point the point in question
     * @param arcLength at which the POCA occurs
     * @param voca vector of closest approach
     * @param docaError error on distance of closest approach
     **/ 
    void pointPocaError(const AcdRecon::TrackData& track, const HepPoint3D& point,
			const double& arcLength, const HepVector3D& voca, double& docaError );
    
    
    
    /**
     * @brief POINT where a track crosses a plane defined by a point and an orientation matrix
     *
     * @param track the track projection data
     * @param planePoint reference point in plane
     * @param norm is the orientation of the plane (normal vector)
     * @param arcLength at which the intersection occurs
     * @param isec is the intersection point
     */
    bool crossesPlane(const AcdRecon::TrackData& track,  
		      const HepPoint3D& planePoint, 
		      const HepVector3D& norm,
		      double& arcLength, HepPoint3D& isec );
    
    /**
     * @brief error on POINT where a track crosses a plane defined by a point and an orientation matrix
     *
     * @param track the track projection data
     * @param planePoint reference point in plane
     * @param toGlobal is the orientation of the plane
     * @param arcLength at which the intersection occurs
     * @param cov projected covarience of intersection
     */
    void crossesPlaneError(const AcdRecon::TrackData& track,  
			   const HepPoint3D& planePoint, 
			   const HepGeom::Transform3D& toLocal,
			   const double& arcLength, 
			   CLHEP::HepSymMatrix& cov );
    
    /**
     * @brief POCA Between a track and a ray defined by two points
     *
     * @param track the track projection data
     * @param c1 first corner
     * @param c2 second corner
     * @param rayDoca data about the POCA
     * @param edgeLen length along the ray between the points to the POCA
     */
    bool rayDoca_withCorner(const AcdRecon::TrackData& track, 
			    const HepPoint3D& c1, const HepPoint3D& c2,
			    double& arcLength, double& rayLength, HepPoint3D& poca, HepVector3D& voca, 
			    int& edgeOrCorner, bool allowCorner = false );
    
    /**
     * @brief Error on POCA Between a track and a ray defined by two points
     *
     * @param track the track projection data
     * @param c1 first corner
     * @param c2 second corner
     * @param rayDoca data about the POCA
     * @param edgeLen length along the ray between the points to the POCA
     */
    void rayDocaError(const AcdRecon::TrackData& track, 
		      const HepPoint3D& c1, const HepPoint3D& c2,
		      const double& arcLength, const double& rayLength, 
		      const HepVector3D& voca, 
		      double& docaErr );
    
    
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
		      AcdRecon::PocaData& data, bool isInside);
    
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
		    AcdRecon::PocaData& data);
    
    
    /**
     * @brief initiliaze the kalman propagator
     *
     * @param prop the propagator
     * @param aTrack the full TKR track representation
     * @param trackData the minimal ACD track representation
     * @param maxArcLength furthers arcLength we want to extend the track to
     */
    //void startPropagator(IPropagator& prop, const Event::TkrTrack& aTrack, const AcdRecon::TrackData& trackData,
    //		 const double& maxArcLength);

    /**
     * @brief initiliaze the kalman propagator
     *
     * @param prop the propagator
     * @param trackParams the track parameters
     * @param trackData the minimal ACD track representation
     * @param maxArcLength furthers arcLength we want to extend the track to
     */
    void startPropagator(IPropagator& prop, const Event::TkrTrackParams& trackParams, const AcdRecon::TrackData& trackData,
			 const double& maxArcLength);
    
    /**
     * @brief run the propagtor out to a specified arclength
     *
     * @param prop the propagator
     * @param arcLength current arclength of the track
     * @param next_params full TKR rep of parameters at next step in Kalman fit
     * @param trackDatat the minimal ACD track representation
     **/
    void propagateToArcLength(IPropagator& prop,
			      const double& arcLength,
			      const AcdRecon::TrackData& trackData,
			      Event::TkrTrackParams& next_params );
    
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

  }
}

#endif
