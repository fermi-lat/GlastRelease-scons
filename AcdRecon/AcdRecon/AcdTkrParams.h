#ifndef ACDTKRPARAMS_H
#define ACDTKRPARAMS_H


#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "Event/Recon/TkrRecon/TkrTrackParams.h"

class AcdTkrParams {
  
  // The track parameters at a point
  
public:
  
  static bool convertToAcdRep(const Event::TkrTrackParams& params, const double& zRef, bool up,
			      HepPoint3D& refPoint, HepVector3D& dir, HepSymMatrix& cov);

public:

  /// C'tor
  AcdTkrParams();
  
  AcdTkrParams(const HepPoint3D& refPoint, const HepVector3D& dir, 
	       const HepSymMatrix& cov);

  AcdTkrParams(const Event::TkrTrackParams& params, const double& zRef, bool up);

  void set(const Event::TkrTrackParams& params, const double& zRef, bool up);

  // Make it look like a vector (params), or like a matrix (covarience), index starts at 1, b/c CLHEP is antiquated
  double operator()(int i) const;
  const double& operator()(int i, int j) const;

  // access
  inline const HepPoint3D& refPoint() const { return m_refPoint; }
  inline const HepVector3D& dir() const { return m_dir; }
  inline const HepSymMatrix& cov() const { return m_cov; }
  
private:
  
  HepPoint3D   m_refPoint;         // the reference point of the track
  HepVector3D  m_dir;              // the direction of the track
  
  HepSymMatrix m_cov;              // 3x3 symetric covarience matrix on the position defined @ the reference point
  
};

#endif
