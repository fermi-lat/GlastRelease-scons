#include "../AcdRecon/AcdTkrParams.h"

#include <stdexcept>

using namespace CLHEP;

bool AcdTkrParams::convertToAcdRep(const Event::TkrTrackParams& params, const double& zRef, bool up,
				   HepPoint3D& refPoint, HepVector3D& dir, HepSymMatrix& cov) {

  // first get the parameters
  // Note that Event::TkrTrackParams does indexing from 1.  I think that this is b/c CLHEP hates us.
  
  // get the reference point
  refPoint.set( params(1), params(2), zRef );

  // get the direction
  // start with the slopes.  Aka tangent.
  const double Mx = params(3);
  const double My = params(4);
  const double tanSqThetaX = Mx*Mx;
  const double tanSqThetaY = My*My;

  // use the relation sin**2 = tan**2 / (1 + tan**2) to get square of directional cosines
  // (directional cosine is sin(theta) = cos(pi/2 - theta) 
  const double sinSqThetaX = tanSqThetaX/ ( 1. + tanSqThetaX );
  const double sinSqThetaY = tanSqThetaY/ ( 1. + tanSqThetaY );
  
  // normalize direction vector to unit length
  const double sinSqThetaZ = 1. - sinSqThetaX - sinSqThetaY;

  // take some square roots
  double xDir = sinSqThetaX > 0 ? sqrt(sinSqThetaX) : 0.;  
  double yDir = sinSqThetaY > 0 ? sqrt(sinSqThetaY) : 0.;
  double zDir = sinSqThetaZ > 0 ? sqrt(sinSqThetaZ) : 0.;
  
  // if this is downgoing flip the direction
  if ( !up) {
    xDir *= -1.;
    yDir *= -1.;
    zDir *= -1.;
  }
  // and set the direction
  dir.set(xDir,yDir,zDir);

  // Now on to the covarience matrix on the position
  // We get this by calculating the intersection of the track 
  // with a plane perpindicular to the track
  
  // the Jacobian is 
  //     1       ,   m_y/ m_x 
  //  m_x/m_y    ,      1
  //    m_x      ,     m_y
  

  double Vxx = cov(1,1);
  double Vxy = cov(1,2);
  double Vyy = cov(2,2);
    
  double one_over_Mx = fabs(Mx) > 1e-6 ? 1./ Mx : 1e-6;
  double one_over_My = fabs(My) > 1e-6 ? 1./ My : 1e-6;
  double Mx2 = Mx*Mx;
  double My2 = My*My;
  double one_over_Mx2 = one_over_Mx * one_over_Mx;
  double one_over_My2 = one_over_My * one_over_My;
  
  cov(1,1) = Vxx + Vxy * (2. * My * one_over_Mx) + Vyy * (My2 * one_over_Mx2);  
  cov(1,2) = Vxx *(Mx * one_over_My) + Vxy * (2.) + Vyy * (My * one_over_Mx);
  cov(1,3) = Vxx *(Mx) + Vxy * (2. * My ) + Vyy * (My2 * one_over_Mx);  
  cov(2,2) = Vxx *(Mx2 * one_over_My2) + Vxy * (2. * Mx * one_over_My) + Vyy;
  cov(2,3) = Vxx *(Mx2 * one_over_My) + Vxy * (2. * Mx ) + Vyy * (My); 
  cov(3,3) = Vxx *(Mx2) + Vxy * (2. * Mx * My) + Vyy * (My2);
  return true;
}

AcdTkrParams::AcdTkrParams()
  :m_refPoint(),    // origin
   m_dir(0.,0.,1.), // vertical track
   m_cov(3,1){      // 3x3 identity matrix  
}

/// C'tor
AcdTkrParams::AcdTkrParams(const HepPoint3D& refPoint, const HepVector3D& dir, 
			   const HepSymMatrix& cov)
  :m_refPoint(refPoint),
   m_dir(dir),
   m_cov(cov){
}

AcdTkrParams::AcdTkrParams(const Event::TkrTrackParams& params, const double& zRef, bool up)
  :m_refPoint(),    // origin
   m_dir(0.,0.,1.), // vertical track
   m_cov(3,1){      // 3x3 identity matrix   
  convertToAcdRep(params,zRef,up,m_refPoint,m_dir,m_cov);
}

void AcdTkrParams::set(const Event::TkrTrackParams& params, const double& zRef, bool up) {
  convertToAcdRep(params,zRef,up,m_refPoint,m_dir,m_cov);  
}

double AcdTkrParams::operator()(int i) const {
  switch (i) {
  case 0:  return m_refPoint.x();
  case 1:  return m_refPoint.y();
  case 2:  return m_refPoint.z();
  case 3:  return m_dir.y();
  case 4:  return m_dir.x();
  case 5:  return m_dir.z();
  }
  throw std::invalid_argument("Invalid index for AcdTkrParams");
  static const double nullDummy(0.);
  return nullDummy;
}

const double& AcdTkrParams::operator()(int i, int j) const {
  return m_cov(i,j);
}
