#ifndef ACDFRAMEUTIL_H
#define ACDFRAMEUTIL_H

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <vector>

/**
 * @brief Some definitions and functions used in handling ACD coordinate frames
 *
 * We have to keep track of 3 frames for each ACD detector element
 * 
 * 1) The global frame, ie, the LAT based frame.  
 *     X, Y run across the towers, Z runs up the towers. 
 *     X, Y = 0 at the center of the LAT, Z = 0 at the TKR-CAL interface
 *
 * 2) The local frame, ie, tile or ribbon frame.
 *     X, Y run across tiles. Z runs out of the tile.
 *     X, Y = 0 at the center of the tile.
 *     for ribbons x runs across ribbon, y runs along ribbon, 
 *                 z out of ribbon plane
 *     X, Y = 0 at the center of the ribbon
 *
 * 3) The intermediate or GEANT frame.
 *    This is the same as the local frame except for any possible tilting of the 
 *    element
 * 
 */

namespace AcdFrameUtil {

  /**
   * @brief Define the various rotations used to describe the ACD
   */

  typedef enum {
    FRAME_NONE = 100000,   // can't find anything suitable
    FRAME_FACE0 = 0,       // natural for top tiles, = LAT global frame
    FRAME_TOP = 0,
    FRAME_XMEAS = 0,       // also natural for x-meas ribbons on top
    FRAME_FACE1 = 1,       // natural for face1 (-X face) tiles
    FRAME_MINUSX =1,
    FRAME_FACE2 = 2,       // natural for face2 (-Y face) tiles
    FRAME_MINUSY = 2,
    FRAME_FACE3 = 3,       // natural for face3 (+X face) tiles
    FRAME_PLUSX = 3,
    FRAME_FACE4 = 4,       // natural for face4 (+Y face) tiles
    FRAME_PLUSY = 4,
    FRAME_FACE1_YDWN = 5,  // like FACE1 but local +Y is global -Z
    FRAME_MINUSX_YDWN = 5,
    FRAME_FACE2_YDWN = 6,  // like FACE2 but local +Y is global -Z
    FRAME_MINUSY_YDWN = 6,
    FRAME_FACE3_YDWN = 7,  // like FACE3 but local +Y is global -Z
    FRAME_PLUSX_YDWN = 7,
    FRAME_FACE4_YDWN = 8,  // like FACE4 but local +Y is global -Z
    FRAME_PLUSY_YDWN = 8,
    FRAME_YMEAS = 9,       // natural for top segments of y-meas ribbons
    FRAME_YMEAS_ZROT180 = 10,  // FRAME_YMEAS rotated 180 deg. about Z
    FRAME_XMEAS_ZROT180 = 11   // FRAME_XMEAS rotated 180 deg. about Z
  } AcdReferenceFrame;
  
  /// pi/2 to lots of accuracy
  const double PI_OVER_2(1.57079632679489656);

  /// No Rotation at all
  const HepGeom::Transform3D s_identity;

  /// Rotations about X 
  const HepGeom::RotateX3D s_rx_pi_over_2(PI_OVER_2);
  const HepGeom::RotateX3D s_rx_pi(2*PI_OVER_2);
  const HepGeom::RotateX3D s_rx_3pi_over_2(3*PI_OVER_2);
  
  /// Rotations about Y
  const HepGeom::RotateY3D s_ry_pi_over_2(PI_OVER_2);
  const HepGeom::RotateY3D s_ry_pi(2*PI_OVER_2);
  const HepGeom::RotateY3D s_ry_3pi_over_2(3*PI_OVER_2);

  /// Rotations about Z
  const HepGeom::RotateZ3D s_rz_pi_over_2(PI_OVER_2);
  const HepGeom::RotateZ3D s_rz_pi(2*PI_OVER_2);
  const HepGeom::RotateZ3D s_rz_3pi_over_2(3*PI_OVER_2);
  
  /**
   * @return the rotation from the GEANT frame to the LOCAL frame
   *
   * This is just a matter of renaming the axes
   */
  const HepGeom::Transform3D& getRotationToLocal(AcdReferenceFrame type);  

  /**
   * @return the rotation from the LOCAL frame to the GEANT frame
   *
   * This is just a matter of renaming the axes
   */
  const HepGeom::Transform3D& getRotationToGeant(AcdReferenceFrame type);

  /**
   * @brief get the the corners in GLOBAL frame
   *
   * @param xVector vector along local X axis, expressed in GLOBAL frame
   * @param yVector vector along local Y axis, expressed in GLOBAL frame
   * @param corner is filled with the locations of the 4 corners.
   *
   * Order is (-,-),(-,+),(+,+),(+,-) 
   *  a.k.a.   clockwise starting from (-,-) corner 
   *
   */  
  void getCornersSquare(const HepPoint3D &center, const HepVector3D& xVector, const HepVector3D& yVector,
			HepPoint3D *corner);
  /**
   * @brief get the the corners in GLOBAL frame
   *
   * @param xVector vector along local X axis at top of tile, expressed in GLOBAL frame
   * @param xVector vector along local X axis at bottom of tile, expressed in GLOBAL frame
   * @param yVector vector along local Y axis, expressed in GLOBAL frame
   * @param corner is filled with the locations of the 4 corners.
   *
   * Order is (-,-),(-,+),(+,+),(+,-) 
   *  a.k.a.   clockwise starting from (-,-) corner 
   *
   */  
  void getCornersTrap(const HepPoint3D &center, const HepVector3D& x1Vector, const HepVector3D& x2Vector,
		      const HepVector3D& yVector,  HepPoint3D *corner);

  // Convert the dimension vector from the way it is expressed in GEANT frame
  // into the LOCAL frame

  /**
   * @brief Convert the dimension vector from the way it is expressed in GEANT frame
   *        into the LOCAL frame
   *
   * @param inGeant dimension vector in GEANT frame
   * @param inLocal dimension vector in LOCAL frame
   */    
  void transformDimensionVector(AcdReferenceFrame type, 
				const std::vector<double>& inGeant, 
				std::vector<double>& inLocal);
  /**
   * @brief Convert the dimension vector from the way it is expressed in GEANT frame
   *        into the LOCAL frame
   *
   * @param inGeant dimension vector in GEANT frame
   * @param inLocal dimension vector in LOCAL frame
   */    
  void transformDimensionVectorTrap(AcdReferenceFrame type, 
				    const std::vector<double>& inGeant, 
				    std::vector<double>& inLocal);


  /**
   * @brief Get the midpoint between two points
   *
   * @param p1 first point
   * @param p2 second point
   * @param mid filled with value of midpoint   
   */    
  void getMidpoint(const HepPoint3D &p1, const HepPoint3D& p2, HepPoint3D& mid);

  /**
   * @brief Build a transform from GEANT to LOCAL frame
   *
   * @param type enum giving orientation of axes
   * @param center of the detector element
   * @param out filled with resulting transform
   */      
  void buildTransform(AcdReferenceFrame type, const HepPoint3D& center,
		      HepGeom::Transform3D& out);
  
  /**
   * @brief Build the composite rotations.  Called once in initialization
   *
   */
  void buildRotations();


  /**
   *
   *
   **/
  void getErrorAxes(const HepGeom::Transform3D& toGlobal, const HepSymMatrix& cov, HepVector3D& v1, HepVector3D& v2);

}

#endif
