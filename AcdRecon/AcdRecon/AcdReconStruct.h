#ifndef ACDRECONSTRUCT_H
#define ACDRECONSTRUCT_H

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "idents/AcdId.h"
#include "idents/AcdGapId.h"
#include "idents/VolumeIdentifier.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "../AcdRecon/AcdGap.h"

#include <map>
#include <list>

namespace AcdRecon {
  
  /// This struct stores the data about the point of closest approach to active elements that we pass around internally
  struct PocaData {
  public:
    PocaData(){reset(2000.);}      
    void reset(double maxDist) {
      m_id = idents::AcdId();
 
      m_arcLengthPlane = 0.;
      m_hitsPlane = Point();
      m_inPlane = Point();
      m_planeError = HepSymMatrix(2,1);
      m_volume = -1;

      m_cosTheta = 0.;
      m_path = 0;
      m_activeX = -maxDist;
      m_activeY = -maxDist;
      m_active2D = -maxDist;

      m_arcLength = 0.;
      m_ribbonLength = 0.;
      m_poca = Point(); 
      m_pocaVector = Vector();
      m_active3D = -maxDist;
      m_active3DErr = 0;
      m_region = -1;

    }
    idents::AcdId m_id;        // The AcdId of the hit element

    double m_arcLengthPlane;   // Length along the track to the plane of the detector
    Point m_hitsPlane;         // 3D point that track crosses detector plane, in global coordiantes
    Point m_inPlane;           // 3D point that track crosses detector plane, in global coordiantes
    HepSymMatrix m_planeError; // 2x2 covarience martix in the detector plane    
    int m_volume;              // Which volume got hit

    double m_cosTheta;        // angle between track and plane normal
    double m_path;            // pathlength of track in plane
    double m_activeX;         // The distance of closest aproach to the relevent X edge in 2D
    double m_activeY;         // The distance of closest aproach to the relevent Y edge in 2D (also, length along ribbon)
    double m_active2D;        // The distance of closest aproach to the relevent edge in 2D

    double m_arcLength;       // Length along the track to the poca
    double m_ribbonLength;    // Length along the ribbon to the poca
    Point m_poca;             // Point of closest approach
    Vector m_pocaVector;      // Vector from Track to POCA

    double m_active3D;        // The distance of closest aproach to the relevent edge in 3D
    double m_active3DErr;     // The error on distance of closest aproach to the relevent edge in 3D

    int m_region;             // One of the enums defined in Event/Recon/AcdRecon/AcdTkrPoca.h 
  };  
  typedef std::map<idents::AcdId,PocaData>         PocaDataMap;
  typedef std::map<idents::AcdId,PocaData*>        PocaDataPtrMap;
 
  /// This struct stores the data about the intersection with the nominal ACD
  struct ExitData {
  public:
    ExitData ()
      :m_face(-1),m_arcLength(-1.){;}
    void reset() {
      m_arcLength = -1.;
      m_face = -1;
    }
    int    m_face;        // 0:top 1:-X 2:-Y 3:+X 4:+Y 5:bottom
    double m_arcLength;   // Length along the track to the m_x
    Point  m_x;           // Intersection Point
  };

  /// this strtuc stores the ACD fiducial volume,  (Should be filled from XML geom svc)
  struct AcdVolume {
  public:
    AcdVolume():
      m_top(754.6),m_sides(840.14),m_bottom(-50.){;}
    AcdVolume(double top, double sides, double bottom):
      m_top(top),m_sides(sides),m_bottom(bottom){;}
    double     m_top;     // top is defined by planes at + 754.6 -> up to stacking of tiles
    double     m_sides;   // sides are defined by planes at +-840.14
    double     m_bottom;  // bottom of the ACD is at the z=-50 plane
  };

  /// this struct stores the track pointing data that we pass around
  struct TrackData {
  public:
    TrackData():m_point(),m_dir(),m_energy(0.),m_index(-1),m_upward(true){;}
    TrackData(const HepPoint3D& point, const HepVector3D& dir, double energy, int index, bool up)
      :m_point(point),m_dir(dir),m_energy(energy),m_index(index),m_upward(up){;}
    HepPoint3D  m_point;       // the start (or end) point of the track
    HepVector3D m_dir;         // the direction of the track
    double      m_energy;      // the energy of the track at the start point
    int         m_index;       // the index number of this track
    bool        m_upward;      // which side of track
  };

  // this is to define if a channel has a hit or not
  typedef std::map<idents::AcdId,unsigned int> AcdHitMap;

  /// This struct stores the data about the backsplash
  struct SplashData {
  public:
    SplashData(){reset(2000.);}      
    void reset(double maxDist) {
      m_id = idents::AcdId();

      m_trackIndex = -1;
      m_region = -1;
      m_calEntryPoint = Point();
      m_calEntryVector = Vector();

      m_tileSolidAngle = -maxDist;
      m_weightedTrackAngle = -maxDist;
      m_weightedPathlength = -maxDist;
    }

    void resetTileData(double maxDist) {
      m_id = idents::AcdId();
      m_tileSolidAngle = -maxDist;
      m_weightedTrackAngle = -maxDist;
      m_weightedPathlength = -maxDist;
    }

    idents::AcdId m_id;       // The AcdId of the hit element   

    /// The index of the associated track
    int m_trackIndex;
    
    /// the region that the track exits from (0 -> cal)
    int m_region;

    /// The Vector from the point the track enters the calorimeter to the tile center
    Point m_calEntryPoint;    
    /// The track vector at the point the track enters the calorimeter
    Vector m_calEntryVector;    
    /// The total solid angle of the tile, seen from the track entry point
    double m_tileSolidAngle;    
    /// The average of the angle between the reconstructed track and the line         
    /// connecting the track entry point in the calorimeter at that point
    /// (weighted by the solid angle of the element)                         
    double m_weightedTrackAngle;    
    /// the average of the pathlength inside the tile along the path from the 
    /// track entry point in the calorimeter to the tile
    /// (weighted by the solid angle of the element)     
    double m_weightedPathlength;
  };  
}

#endif
