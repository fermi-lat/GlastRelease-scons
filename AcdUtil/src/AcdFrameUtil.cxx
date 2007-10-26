

#include "../AcdUtil/AcdFrameUtil.h"

namespace AcdFrameUtil {
    
  //const double PI_OVER_2 = 1.57079632679489656;

  std::vector<HepGeom::Transform3D> _toLocal;
  std::vector<HepGeom::Transform3D> _toGeant;

  const HepGeom::Transform3D& getRotationToLocal(AcdReferenceFrame type) {
    // return the rotation from the GLOBAL frame to the LOCAL frame
    // This is just a matter of renaming the axes
    if ( type == FRAME_NONE ) { return _identity; }
    if ( _toLocal.size() == 0 ) {
      buildRotations();
    }
    return _toLocal[type];
  }

  const HepGeom::Transform3D& getRotationToGeant(AcdReferenceFrame type) {
    // return the rotation from the LOCAL frame to the GLOBAL frame
    // This is just a matter of renaming the axes 
    if ( type == FRAME_NONE ) { return _identity; }
    if ( _toGeant.size() == 0 ) {
      buildRotations();
    }
    return _toGeant[type];
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
    corner[0] = center - x2Vector - yVector;
    corner[1] = center - x1Vector + yVector;
    corner[2] = center + x1Vector + yVector;
    corner[3] = center + x2Vector - yVector;
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
    // Order is X1, X2, X_diff, Y, Z
    inLocal[0] = (r.xx() * inGeant[0]) + (r.xy() * inGeant[3]) + (r.xz() * inGeant[3]);
    inLocal[1] = (r.xx() * inGeant[1]) + (r.xy() * inGeant[3]) + (r.xz() * inGeant[3]);
    inLocal[2] = (r.xx() * inGeant[2]) + (r.xy() * inGeant[3]) + (r.xz() * inGeant[3]);
    inLocal[3] = (r.yx() * inGeant[0]) + (r.yy() * inGeant[3]) + (r.yz() * inGeant[3]);
    inLocal[4] = (r.zx() * inGeant[0]) + (r.zy() * inGeant[3]) + (r.zz() * inGeant[3]);
  
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

    _toLocal.clear();
    _toLocal.resize(12);
    _toGeant.clear();
    _toGeant.resize(12);

    // 0 is identity
    //_toLocal[0];
    //_toGeant[0];
    
    // 1 is -X face
    _toLocal[1] = _ry_pi_over_2* _rx_3pi_over_2;    
    _toGeant[1] = _rx_pi_over_2* _ry_3pi_over_2;     
    
    // 2 is -Y face
    _toLocal[2] = _rx_3pi_over_2;    
    _toGeant[2] = _rx_pi_over_2;    

    // 3 is +X face
    _toLocal[3] = _ry_3pi_over_2*_rx_3pi_over_2;    
    _toGeant[3] = _rx_pi_over_2*_ry_pi_over_2;

    // 4 is +Y face
    _toLocal[4] = _ry_pi*_rx_3pi_over_2;  
    _toGeant[4] = _rx_pi_over_2*_ry_pi;

    // 5 is -X face, rotated by 180 so that y goes down
    _toLocal[5] = _rz_pi*_toLocal[1];
    _toGeant[5] = _toGeant[1]*_rz_pi;
   
    // 5 is -Y face, rotated by 180 so that y goes down
    _toLocal[6] = _rz_pi*_toLocal[2];
    _toGeant[6] = _toGeant[2]*_rz_pi;
    
    // 7 is -X face, rotated by 180 so that y goes down
    _toLocal[7] = _rz_pi*_toLocal[3];
    _toGeant[7] = _toGeant[3]*_rz_pi;
    
    // 8 is -Y face, rotated by 180 so that y goes down
    _toLocal[8] = _rz_pi*_toLocal[4];
    _toGeant[8] = _toGeant[4]*_rz_pi;
        
    // 9 is top. rotated by 90 for Y-ribbons
    _toLocal[9] = _rz_pi_over_2;    
    _toGeant[9] = _rz_3pi_over_2;

    // 10 is top. rotated by -90 for Y-ribbons in opposite direction (not used)
    _toLocal[10] = _rz_3pi_over_2;
    _toGeant[10] = _rz_pi_over_2;
   
    // 11 is top. rotated by 180 for X-ribbons in opposite direction (not used)
    _toLocal[11] = _rz_pi;
    _toGeant[11] = _rz_pi;
    
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
}
