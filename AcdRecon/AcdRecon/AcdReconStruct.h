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

namespace Event {
  class AcdHit;
}

namespace AcdRecon {

  /**
   *  @struct PocaData
   *
   *  @brief This struct stores the data about the point of closest 
   *  approach to active elements that we pass around internally
   *
   */
  struct PocaData {
  public:
    /// c'tor just nulls data
    PocaData(){reset(2000.);} 

    /// Reset all the value to defaults (either 0. or -2000)
    void reset(double maxDist) {
      m_id = idents::AcdId();
 
      m_arcLengthPlane = 0.;
      m_hitsPlane = HepPoint3D();
      m_inPlane = HepPoint3D();
      m_planeError_proj = HepSymMatrix(2,0);
      m_planeError_prop = HepSymMatrix(2,0);
      m_volumePlane = -1;

      m_cosTheta = 0.;
      m_path = 0;
      m_activeX = -maxDist;
      m_activeY = -maxDist;
      m_active2D = -maxDist;

      m_arcLength = 0.;
      m_rayLength = 0.;
      m_ribbonLength = 0.;
      m_poca = HepPoint3D(); 
      m_voca = HepVector3D();
      m_active3D = -maxDist;
      m_active3DErr_proj = 0;
      m_active3DErr_prop = 0;
      m_volume = -1;
      m_region = -1;

      m_hasHit = false;
    }
    /// The AcdId of the hit element
    idents::AcdId m_id;       

    /// Length along the track to the plane of the detector
    double m_arcLengthPlane;   
    /// 3D point that track crosses detector plane, in global coordiantes
    HepPoint3D m_hitsPlane;    
    /// 3D point that track crosses detector plane, in local coordiantes
    HepPoint3D m_inPlane; 
    /// 2x2 covarience martix in the detector plane, projected from original track params
    HepSymMatrix m_planeError_proj; 
    /// 2x2 covarience martix in the detector plane, using propagated error
    HepSymMatrix m_planeError_prop;     
    /// Which volume got hit, 
    int m_volumePlane;             

    /// angle between track and plane normal
    double m_cosTheta;      
    /// pathlength of track in plane
    double m_path;    
    /// The distance of closest aproach to the relevent X edge in 2D        
    double m_activeX;    
    /// The distance of closest aproach to the relevent Y edge in 2D (also, length along ribbon)      
    double m_activeY;
    /// The distance of closest aproach to the relevent edge in 2D  
    double m_active2D;        

    /// Length along the track to the poca
    double m_arcLength;      
    /// Length along the ribbon segment or tile edge to the poca
    double m_rayLength;  
    /// Length along the full ribbon
    double m_ribbonLength;
    /// Point of closest approach
    HepPoint3D m_poca;
    /// Vector from Track to POCA
    HepVector3D m_voca;

    /// The distance of closest aproach to the relevent edge in 3D
    double m_active3D;        
    /// The error on distance of closest aproach to the relevent edge in 3D, projected from original track params
    double m_active3DErr_proj;     
    /// The error on distance of closest aproach to the relevent edge in 3D, using propagated error
    double m_active3DErr_prop;     

    /// Which volue was used to get poca
    int m_volume;
    /// Which edge or corner was used to get poca 
    int m_region;

    /// Does the poca have a hit
    bool m_hasHit;
  };  

  /// Define maps from AcdId to PocaData
  typedef std::map<idents::AcdId,PocaData>         PocaDataMap;
  typedef std::map<idents::AcdId,PocaData*>        PocaDataPtrMap;
 
  /**
   *  @struct ExitData
   *
   *  @brief This struct stores the data about the intersection with the nominal ACD
   *  (Basically a box with the dimensions of the ACD, without all the details about 
   *  shingling and gaps)
   *  
   */
  struct ExitData {
  public:
    /// c'tor just nulls data
    ExitData ()
      :m_face(-1),m_arcLength(-1.),m_arcTol(0.){;}

    /// Reset all the value to defaults (-1)
    void reset() {
      m_arcLength = -1.;
      m_arcTol = 0.;
      m_face = -1;
    }
    /// which face of ACD 0:top 1:-X 2:-Y 3:+X 4:+Y 5:bottom
    int    m_face;        
    /// Length along the track to the m_x
    double m_arcLength;   
    /// Possible error in pathlength
    double m_arcTol;
    /// Intersection Point
    Point  m_x;           

    /// 2x2 covarience martix in the ACD plane, projected from original track params
    HepSymMatrix m_planeError_proj; 
    /// 2x2 covarience martix in the ACD plane, using propagated error
    HepSymMatrix m_planeError_prop;     

  };

  /**
   *  @struct AcdVolume
   *
   *  @brief this struct store the information about the nominal size of the ACD
   *  
   */
  struct AcdVolume {
  public:

    /// null c'tor give default values
    AcdVolume()
      :m_top(754.6),m_sides(840.14),m_bottom(-50.),
       m_topTol(30.),m_sideTol(30.),m_botTol(0.){;}
    
    /// standard c'tor allows setting values
    AcdVolume(double top, double sides, double bottom,
	      double topTol, double sideTol, double botTol)
      :m_top(top),m_sides(sides),m_bottom(bottom),
       m_topTol(topTol),m_sideTol(sideTol),m_botTol(botTol)   {;}
    
    /// top is defined by planes at + 754.6 -> up to stacking of tiles
    double     m_top;     
    /// sides are defined by planes at +-840.14
    double     m_sides;   
    /// bottom of the ACD is at the z=-50 plane
    double     m_bottom;  
    /// tolerence at top +- 20mm 
    double     m_topTol;     
    /// tolerence at sides +-15mm
    double     m_sideTol;   
    /// tolerence at bottom is 0
    double     m_botTol;  
  };

  /**
   *  @struct TrackData
   *
   *  @brief this struct stores the track pointing data that we pass around
   *  
   */
  struct TrackData {
  public:
    /// Null c'tor, just a place keeper
    TrackData():
      m_point(),m_dir(),m_energy(0.),m_index(-1),m_upward(true),
      m_cov_orig(5,0),m_current(),m_cov_prop(5,0){;}
    /// C'tor from point, direction, energy, track index, up or down flag
    TrackData(const HepPoint3D& point, const HepVector3D& dir, 
	      double energy, int index, bool up)
      :m_point(point),m_dir(dir),m_energy(energy),m_index(index),m_upward(up),
       m_cov_orig(5,0),m_current(point),m_cov_prop(5,0){;}
    
    /// the start (or end) point of the track
    HepPoint3D  m_point;       
    /// the direction of the track
    HepVector3D m_dir;         
    /// the energy of the track at the start point (needed for Kalman propagator)
    double      m_energy;      
    /// the index number of this track
    int         m_index;       
    /// which side of track
    bool        m_upward;  

    /// 5x5 covarience martix in the ACD rep, from original track params
    HepSymMatrix m_cov_orig; 

    /// The projected point of the track
    mutable HepPoint3D  m_current;  

    /// 5x5 covarience martix in the ACD rep, propagated error
    mutable HepSymMatrix m_cov_prop;     

  };

  /// this is to define if a channel has a hit or not
  typedef std::map<idents::AcdId,const Event::AcdHit*> AcdHitMap;


  /**
   *  @struct SplashData
   *
   *  @brief  This struct stores the data about the backsplash
   *  
   */
  struct SplashData {
  public:

    /// C'tor just nulls all the data
    SplashData(){reset(2000.);}  

    /// Reset all the value to defaults (either -1. or -2000)
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

    /// Reset only the tile related data to defaults (-2000)
    void resetTileData(double maxDist) {
      m_id = idents::AcdId();
      m_tileSolidAngle = -maxDist;
      m_weightedTrackAngle = -maxDist;
      m_weightedPathlength = -maxDist;
    }

    /// The AcdId of the hit element   
    idents::AcdId m_id;       

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
