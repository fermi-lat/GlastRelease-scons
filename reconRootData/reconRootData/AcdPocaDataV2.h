#ifndef ROOT_ACDPOCADATAV2_H
#define ROOT_ACDPOCADATAV2_H

#include "TObject.h"
/** 
 * @class AcdPocaDataV2
 *
 * @brief ROOT object information about the Point of Closest Approach (POCA) between an extrapolated track
 *  and the Acd (either a detector element or a gap).  
 *
 *  This POCA is calculated in 3D.  The doca is defined to be positive if the track goes inside the active distance 
 *  or gap and negative otherwise.  Both AcdTkrHitPoca and AcdTkrGapPoca inherit from this class.
 *
 *  The main access functions are:
 *    - float getArcLength()
 *      - which returns the arc-length along the track at which the POCA occurs.  
 *        Postive values are given for intersections above the first track hit.
 *    - float getDoca(), float getDocaErr()
 *      - which returns the distance of closest approach.  This is usually expressed as an active distance
 *        Positive values mean the track went into the volume in question, 
 *        Negative values mean it missed the volume.
 *        This magnitude is the distance to the edge of the volume.
 *    - const TVector3& getPoca()
 *      - which returns the Point of Closest Approach (POCA) in global coordinates
 *    - const TVector3& getPocaVector()
 *      - which returns the  vector between the POCA and the closest edge of the volume in question
 *
 * @author Eric Charles
 *
 * $Header$
 */

#include "TVector3.h"

class AcdPocaDataV2 {
  
public:

  /// Default constructor.  
  AcdPocaDataV2();
  
  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdPocaDataV2(Int_t volume, Int_t region, Float_t arcLength, 
              Float_t doca, Float_t docaErrProj,  Float_t docaErrProp,
              const TVector3& poca, const TVector3& pocaVector);

  /// Copy constructor
  AcdPocaDataV2(const AcdPocaDataV2& other);

  virtual ~AcdPocaDataV2() {;}
  
  /// Assignment operator
  AcdPocaDataV2& operator=(const AcdPocaDataV2& other);
    
  inline Int_t getVolume() const {
    return m_volume;
  }

  inline Int_t getRegion() const {
    return m_region;
  }

  /// Return the arcLength along the track at which the POCA occurs
  /// This is calculated in 3D.  See AcdPocaTool for details  
  inline Float_t getArcLength() const {
    return m_arcLength;
  }

  /// Return the distance of closest approach (doca)
  /// This is calculated in 3D.  See AcdPocaTool for details  
  inline Float_t getDoca() const {
    return m_doca;
  }

  /// Return the error on distance of closest approach (doca)
  /// This is calulated the full Kalman Filter propagation of the track
  /// to the POCA and the projection of the propagated covarience matrix 
  /// along the line of the doca.
  inline Float_t getDocaErrProj() const {
    return m_docaErr_proj;
  }

  /// Return the error on distance of closest approach (doca)
  /// This is calulated the full Kalman Filter propagation of the track
  /// to the POCA and the projection of the propagated covarience matrix 
  /// along the line of the doca.
  inline Float_t getDocaErrProp() const {
    return m_docaErr_prop;
  }

  
  
  /// Return the POCA (in global coordinates)
  const TVector3& getPoca() const {
    return m_poca;
  }
 
  /// Return vector between the POCA and the closest edge of the detector element or GAP 
  const TVector3& getPocaVector() const {
    return m_pocaVector;
  }
 
  
  /// set all the values
  void set(Int_t volume, Int_t region, Float_t arcLength,
           Float_t doca, Float_t docaErrProj, Float_t docaErrProp,
           const TVector3& poca, const TVector3& pocaVector);

  /// copy in the object
  void copy(const AcdPocaDataV2& other);

  /// set individaul values

  // for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdPocaDataV2&, const std::string & name = "" ) const ; // for tests

  /// reset all the values to their default (lowercase to avoid conflict with TObject::Clear())
  void clear();
  
  /// Print out this structure (lowercase to avoid conflict with TObject::Plear())
  void print(Option_t *option) const;
  
private:
  
  /// Which volume of the Tile or Ribbon does this occur w.r..t
  Int_t   m_volume;

  /// Which region of the poca occur w.r.t. (ie, which edge or corner of the tile or ribbon)
  Int_t   m_region;
  
  /// The arclength at which the poca occurs
  Float_t m_arcLength;

  /// The Distance of Closest Approach between track and tile or ribbon
  Float_t m_doca;                 
  
  /// The Error on the DOCA, this is the track error ellipsoid project along the doca direction
  Float_t m_docaErr_proj;

  /// The Error on the DOCA, this is the track error ellipsoid project along the doca direction
  Float_t m_docaErr_prop;  
  
  /// The Point of Closest Approach
  TVector3 m_poca;

  /// The vector from the POCA to the Detector Element edge
  TVector3 m_pocaVector;
  
  ClassDef(AcdPocaDataV2,1)
  
};
#endif
