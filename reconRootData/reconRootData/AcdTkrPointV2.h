#ifndef ROOT_ACDTKRPOINTV2_H
#define ROOT_ACDTKRPOINTV2_H

#include "TObject.h"

#include "AcdTkrLocalCoordsV2.h"

#include "TVector3.h"

/**
*  @class AcdTkrPointV2
*
*  @brief ROOT object which stores information where an extrapolated track exits the nominal ACD volume.  
*
*  This Point is calculated in 3D.  We include information the track paramterization at the point.
*  
*  The main access functions are:
*    - float getArcLength()
*      - which returns the arc-length along the track at which the point occurs.  
*        Postive values are given for intersections above the first track hit.
*    - const TVector3& point()
*      - which returns the intersection point in global coordinates
*    - int face() 
*      - which return the side of the ACD did we exited (0=Top, 1=-X, 2=-Y, 3=+X, 4=+Y, 5=Bottom)
*    - const Event::TkrTrackParams& paramsAtPoint()
*      - which returns the kalman propagated track parameters at the POINT
*  
*  \author Eric Charles
*
* $Header$
*/

  
class AcdTkrPointV2 : public TObject, public AcdTkrLocalCoordsV2 {
  
public:
  
  /// Default constructor.  
  AcdTkrPointV2();
  
  /// Copy constructor
  AcdTkrPointV2(const AcdTkrPointV2& other);

  /// Destructor is trivial
  virtual ~AcdTkrPointV2() {;}

  /// Assignment operator
  AcdTkrPointV2& AcdTkrPointV2::operator=(const AcdTkrPointV2& other);
    
  inline Int_t getTrackIndex() const {
    return m_trackIndex;
  }
  
  /// set all the values at once
  void set(Int_t trackIndex);
  
  // for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdTkrPointV2&, const std::string & name = "" ) const ; // for tests

  /// reset all the values to their default (lowercase to avoid conflict with TObject::Clear())
  void Clear(Option_t *option);
  
  /// Print out this structure (lowercase to avoid conflict with TObject::Clear())
  void Print(Option_t *option) const;

private:
    
  /// Which track
  Int_t    m_trackIndex;
    
  ClassDef(AcdTkrPointV2,1)

};

#endif
