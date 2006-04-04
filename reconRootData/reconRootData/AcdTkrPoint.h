#ifndef ROOT_ACDTKRPOINT_H
#define ROOT_ACDTKRPOINT_H

#include "TObject.h"

#include "TkrTrackParams.h"

#include "TVector3.h"

/**
*  @class AcdTkrPoint
*
*
*  @brief This class stores information about the Point of Closest Approach (POINT) between an extrapolated track
*  and a hit Acd element (tile or ribbon).  This POINT is calculated in 3D.  The doca is defined to be positive 
*  if the track goes inside the active distance and negative otherwise
*  
*  \author Eric Charles
*
* $Header$
*/

  
class AcdTkrPoint : public TObject {
  
public:
  
  /// Default constructor.  
  AcdTkrPoint();
  
  /// Constructor for use in persistent -> transient conversion and reconstruction
  /// Takes arguements as they are stored in ROOT and caluclated by AcdReconAlg
  AcdTkrPoint(double arcLength, int face, 
	      const TVector3& point, const TkrTrackParams& paramsAtPoint);

  /// Copy constructor
  AcdTkrPoint(const AcdTkrPoint& other);

  /// Destructor is trivial
  virtual ~AcdTkrPoint() {;}

  /// Assignment operator
  AcdTkrPoint& AcdTkrPoint::operator=(const AcdTkrPoint& other);
  
  /// Return the arclength along the track at which the track cross the nominal ACD volume
  /// This is calculated in 3D.  
  inline Float_t getArcLength() const {
    return m_arcLength;
  }
  
  /// return the face of the 
  inline int getFace() const {
    return (int)(m_face);
  }
  
  /// Return the POINT (in global coordinates)
  inline const TVector3& getPoint() const {
    return m_point;
  }
  
  /// Return the kalman propagated track parameters at the POINT
  inline const TkrTrackParams& getParamsAtPoint() const {
    return m_paramsAtPoint;
  }
  
  /// set all the values at once
  void set(double arcLength, int face, 
	   const TVector3& point, const TkrTrackParams& paramsAtPoint);
  
  // set the individual values
  inline void setArcLength(float val) { m_arcLength = val; };
  inline void setFace(int val) { m_face = (Char_t)(val); };
  inline void setPoint(const TVector3& val) { m_point = val; };
  inline void setParams(const TkrTrackParams& val) { m_paramsAtPoint = val; };
  
  // for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdTkrPoint&, const std::string & name = "" ) const ; // for tests

  /// reset all the values to their default (lowercase to avoid conflict with TObject::Clear())
  void Clear(Option_t *option);
  
  /// Print out this structure (lowercase to avoid conflict with TObject::Clear())
  void Print(Option_t *option) const;

private:
    
  /// The arclength at which this track crosses the nominal ACD volume
  Float_t m_arcLength;
  
  /// which face of the ACD does the track leaves by
  UChar_t m_face;
  
  /// The point of closest approach of the track to the tile or ribbon
  TVector3 m_point;
  
  /// The track parameters (and covarience matrix) at the POINT
  TkrTrackParams m_paramsAtPoint;

  ClassDef(AcdTkrPoint,1)

};

#endif
