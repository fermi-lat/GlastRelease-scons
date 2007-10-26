#ifndef ACDFRAMEUTIL_H
#define ACDFRAMEUTIL_H

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Geometry/Point3D.h"

#include <vector>


namespace AcdFrameUtil {

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
  
  // pi/2 to lots of accuracy
  const double PI_OVER_2(1.57079632679489656);

  const HepGeom::Transform3D _identity;

  // Rotations about X 
  const HepGeom::RotateX3D _rx_pi_over_2(PI_OVER_2);
  const HepGeom::RotateX3D _rx_pi(2*PI_OVER_2);
  const HepGeom::RotateX3D _rx_3pi_over_2(3*PI_OVER_2);
  
  // Rotations about Y
  const HepGeom::RotateY3D _ry_pi_over_2(PI_OVER_2);
  const HepGeom::RotateY3D _ry_pi(2*PI_OVER_2);
  const HepGeom::RotateY3D _ry_3pi_over_2(3*PI_OVER_2);

  // Rotations about Z
  const HepGeom::RotateZ3D _rz_pi_over_2(PI_OVER_2);
  const HepGeom::RotateZ3D _rz_pi(2*PI_OVER_2);
  const HepGeom::RotateZ3D _rz_3pi_over_2(3*PI_OVER_2);
  
  // return the rotation from the GEANT frame to the LOCAL frame
  // This is just a matter of renaming the axes
  const HepGeom::Transform3D& getRotationToLocal(AcdReferenceFrame type);  

  // return the rotation from the LOCAL frame to the GEANT frame
  // This is just a matter of renaming the axes 
  const HepGeom::Transform3D& getRotationToGeant(AcdReferenceFrame type);

  // Return the corners in GLOBAL frame
  // given vectors along the local X and Y axes, expressed in GLOBAL frame
  // Order is (-,-),(-,+),(+,+),(+,-) 
  // a.k.a.   clockwise starting from (-,-) corner  
  void getCornersSquare(const HepPoint3D &center, const HepVector3D& xVector, const HepVector3D& yVector,
			HepPoint3D *corner);

  // Return the corners in GLOBAL frame
  // given vectors along the local X and Y axes, expressed in GLOBAL frame
  // Order is (-,-),(-,+),(+,+),(+,-) 
  // a.k.a.   clockwise starting from (-,-) corner  
  void getCornersTrap(const HepPoint3D &center, const HepVector3D& x1Vector, const HepVector3D& x2Vector,
		      const HepVector3D& yVector,  HepPoint3D *corner);

  // Convert the dimension vector from the way it is expressed in GEANT frame
  // into the LOCAL frame
  void transformDimensionVector(AcdReferenceFrame type, 
				const std::vector<double>& inGeant, 
				std::vector<double>& inLocal);

  // Convert the dimension vector from the way it is expressed in GEANT frame
  // into the LOCAL frame
  void transformDimensionVectorTrap(AcdReferenceFrame type, 
				    const std::vector<double>& inGeant, 
				    std::vector<double>& inLocal);


  // Get the mid point between two points
  void getMidpoint(const HepPoint3D &p1, const HepPoint3D& p2, HepPoint3D& mid);

  void buildTransform(AcdReferenceFrame type, const HepPoint3D& center,
		      HepGeom::Transform3D& out);
  

  void buildRotations();

}

#endif
