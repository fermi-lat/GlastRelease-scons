

#include "../AcdUtil/AcdFrameUtil.h"
#include "CLHEP/Matrix/Matrix.h"

namespace AcdFrameUtil {
    
  //const double PI_OVER_2 = 1.57079632679489656;

  std::vector<HepGeom::Transform3D>* s_pToLocal = 0;
  std::vector<HepGeom::Transform3D>* s_pToGeant = 0;

  const HepGeom::Transform3D& getRotationToLocal(AcdReferenceFrame type) {
    // return the rotation from the GLOBAL frame to the LOCAL frame
    // This is just a matter of renaming the axes
    if (!s_pToLocal) s_pToLocal = new std::vector<HepGeom::Transform3D>;
    if ( type == FRAME_NONE ) { return s_identity; }
    if ( s_pToLocal->size() == 0 ) {
      buildRotations();
    }
    return (*s_pToLocal)[type];
  }

  const HepGeom::Transform3D& getRotationToGeant(AcdReferenceFrame type) {
    // return the rotation from the LOCAL frame to the GLOBAL frame
    // This is just a matter of renaming the axes 
    if ( type == FRAME_NONE ) { return s_identity; }
    if (!s_pToGeant) s_pToGeant = new std::vector<HepGeom::Transform3D>;
    if ( s_pToGeant->size() == 0 ) {
      buildRotations();
    }
    return (*s_pToGeant)[type];
  }

  void getCornersSquare(const HepPoint3D &center, const HepVector3D& xVector, const HepVector3D& yVector,
			HepPoint3D *corner) {    
    // The corners are returned in order (-,-), (-,+), (+,+), (+,-)
    // a.k.a.   clockwise starting from (-,-) corner
    corner[0] = center - xVector - yVector; 
    corner[1] = center - xVector + yVector;
    corner[2] = center + xVector + yVector;
    corner[3] = center + xVector - yVector;
  }

  void getCornersTrap(const HepPoint3D &center, const HepVector3D& x1Vector, const HepVector3D& x2Vector,
		      const HepVector3D& yVector,  HepPoint3D *corner) {    
    // The corners are returned in order (-,-), (-,+), (+,+), (+,-)
    // a.k.a.   clockwise starting from (-,-) corner
    corner[0] = center - x1Vector - yVector;
    corner[1] = center - x2Vector + yVector;
    corner[2] = center + x2Vector + yVector;
    corner[3] = center + x1Vector - yVector;
  }
  

  void transformDimensionVector(AcdReferenceFrame type,
				const std::vector<double>& inGeant,
				std::vector<double>& inLocal) {
    // Convert the dimension vector from the way it is expressed in GEANT frame
    // into the LOCAL frame
    HepVector3D v(inGeant[0],inGeant[1],inGeant[2]);
    HepVector3D vOut = getRotationToLocal(type)*v;
    inLocal.clear();
    inLocal.push_back(fabs(vOut.x()));
    inLocal.push_back(fabs(vOut.y()));
    inLocal.push_back(fabs(vOut.z()));    
  }

  // Convert the dimension vector from the way it is expressed in GEANT frame
  // into the LOCAL frame
  void transformDimensionVectorTrap(AcdReferenceFrame type, 
				    const std::vector<double>& inGeant, 
				    std::vector<double>& inLocal) {
    // Convert the dimension vector from the way it is expressed in GEANT frame
    // into the LOCAL frame
    const HepGeom::Transform3D& r = getRotationToLocal(type);

    inLocal.resize(5);
    // Input Order is X1, X2, X_diff, Y, Z
    // Output order is X1, Y, Z, X2, X_Diff
    inLocal[0] = (r.xx() * inGeant[0]) + (r.xy() * inGeant[3]) + (r.xz() * inGeant[4]);
    inLocal[1] = (r.yx() * inGeant[0]) + (r.yy() * inGeant[3]) + (r.yz() * inGeant[4]);
    inLocal[2] = (r.zx() * inGeant[0]) + (r.zy() * inGeant[3]) + (r.zz() * inGeant[4]);  
    inLocal[3] = (r.xx() * inGeant[1]) + (r.xy() * inGeant[3]) + (r.xz() * inGeant[4]);
    inLocal[4] = (r.xx() * inGeant[2]) + (r.xy() * inGeant[3]) + (r.xz() * inGeant[4]);
  }

  // Get the mid point between two points
  void getMidpoint(const HepPoint3D &p1, const HepPoint3D& p2, HepPoint3D& mid) {
    HepVector3D v = p1 - p2;
    v *= 0.5;
    mid = p1 - v;
  }


  void buildTransform(AcdReferenceFrame type, const HepPoint3D& center,
		      HepGeom::Transform3D& out){
    
    HepVector3D v(-center.x(),-center.y(),-center.z());
    const HepGeom::Translate3D t(v);
    const HepGeom::Transform3D& r = getRotationToLocal(type);
    out = r*t;
  }
   

  void buildRotations() {

    if (!s_pToLocal) s_pToLocal = new std::vector<HepGeom::Transform3D>;
    if (!s_pToGeant) s_pToGeant = new std::vector<HepGeom::Transform3D>;
    s_pToLocal->clear();
    s_pToLocal->resize(12);
    s_pToGeant->clear();
    s_pToGeant->resize(12);

    // 0 is identity
    //_toLocal[0];
    //_toGeant[0];
    
    // 1 is -X face
    (*s_pToLocal)[1] = s_ry_pi_over_2* s_rx_3pi_over_2;    
    (*s_pToGeant)[1] = s_rx_pi_over_2* s_ry_3pi_over_2;     
    
    // 2 is -Y face
    (*s_pToLocal)[2] = s_rx_3pi_over_2;    
    (*s_pToGeant)[2] = s_rx_pi_over_2;    

    // 3 is +X face
    (*s_pToLocal)[3] = s_ry_3pi_over_2*s_rx_3pi_over_2;    
    (*s_pToGeant)[3] = s_rx_pi_over_2*s_ry_pi_over_2;

    // 4 is +Y face
    (*s_pToLocal)[4] = s_ry_pi*s_rx_3pi_over_2;  
    (*s_pToGeant)[4] = s_rx_pi_over_2*s_ry_pi;

    // 5 is -X face, rotated by 180 so that y goes down
    (*s_pToLocal)[5] = s_rz_pi*(*s_pToLocal)[1];
    (*s_pToGeant)[5] = (*s_pToGeant)[1]*s_rz_pi;
   
    // 5 is -Y face, rotated by 180 so that y goes down
    (*s_pToLocal)[6] = s_rz_pi*(*s_pToLocal)[2];
    (*s_pToGeant)[6] = (*s_pToGeant)[2]*s_rz_pi;
    
    // 7 is -X face, rotated by 180 so that y goes down
    (*s_pToLocal)[7] = s_rz_pi*(*s_pToLocal)[3];
    (*s_pToGeant)[7] = (*s_pToGeant)[3]*s_rz_pi;
    
    // 8 is -Y face, rotated by 180 so that y goes down
    (*s_pToLocal)[8] = s_rz_pi*(*s_pToLocal)[4];
    (*s_pToGeant)[8] = (*s_pToGeant)[4]*s_rz_pi;
        
    // 9 is top. rotated by 90 for Y-ribbons
    (*s_pToLocal)[9] = s_rz_pi_over_2;    
    (*s_pToGeant)[9] = s_rz_3pi_over_2;

    // 10 is top. rotated by -90 for Y-ribbons in opposite direction (not used)
    (*s_pToLocal)[10] = s_rz_3pi_over_2;
    (*s_pToGeant)[10] = s_rz_pi_over_2;
   
    // 11 is top. rotated by 180 for X-ribbons in opposite direction (not used)
    (*s_pToLocal)[11] = s_rz_pi;
    (*s_pToGeant)[11] = s_rz_pi;
    
    /*
    //for ( int i(0); i < 12; i++ ) {
    //std::cout << i << std::endl;
    //std::cout << (int)_toLocal[i](0,0) << ' ' << (int)_toLocal[i](0,1) << ' ' << (int)_toLocal[i](0,2 )<< std::endl;
    //std::cout << (int)_toLocal[i](1,0) << ' ' << (int)_toLocal[i](1,1) << ' ' << (int)_toLocal[i](1,2 )<< std::endl;
    //std::cout << (int)_toLocal[i](2,0) << ' ' << (int)_toLocal[i](2,1) << ' ' << (int)_toLocal[i](2,2 )<< std::endl;
    //std::cout << std::endl;
    //std::cout << (int)_toGeant[i](0,0) << ' ' << (int)_toGeant[i](0,1) << ' ' << (int)_toGeant[i](0,2 )<< std::endl;
    //std::cout << (int)_toGeant[i](1,0) << ' ' << (int)_toGeant[i](1,1) << ' ' << (int)_toGeant[i](1,2 )<< std::endl;
    //std::cout << (int)_toGeant[i](2,0) << ' ' << (int)_toGeant[i](2,1) << ' ' << (int)_toGeant[i](2,2 )<< std::endl;
    //std::cout << std::endl;
    //HepGeom::Transform3D check = _toGeant[i]*_toLocal[i];
    //std::cout << check(0,0) << ' ' << check(0,1) << ' ' << check(0,2 )<< std::endl;
    //std::cout << check(1,0) << ' ' << check(1,1) << ' ' << check(1,2 )<< std::endl;
    //std::cout << check(2,0) << ' ' << check(2,1) << ' ' << check(2,2 )<< std::endl;
    //} 
    */

  }
  
  void getErrorAxes(const HepGeom::Transform3D& toGlobal, const HepSymMatrix& cov, HepVector3D& v1, HepVector3D& v2) {
    
    // First get the rotation axis in the plane
    double top = -2. * cov(1,2);
    double bot = cov(1,1) - cov(2,2);
    double angle = 0.5 * atan2(top,bot);
    
    // Set up the rotation matrix
    HepMatrix rot(2,2,0);
    rot(1,1) = sin(angle);
    rot(2,2) = sin(angle);
    rot(2,1) = cos(angle);
    rot(1,2) = -1. * cos(angle);

    // Rotate out to the get eigenvalues
    HepSymMatrix covUV = cov.similarity(rot);

    // Get the error axis in the UV frame
    HepSymMatrix errUV(2,0);
    errUV(1,1) = sqrt( covUV(1,1) );
    errUV(2,2) = sqrt( covUV(2,2) );

    // Rotate back to the the error in the plane
    HepMatrix errPlane = rot.T() * errUV;
    
    // make the err axis in the plane
    HepVector3D v1l( errPlane(1,1), errPlane(1,2), 0. );
    HepVector3D v2l( errPlane(2,1), errPlane(2,2), 0. );
    
    v1 = toGlobal*v1l;
    v2 = toGlobal*v2l;

  }



}
