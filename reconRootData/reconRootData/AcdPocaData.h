#ifndef ROOT_ACDPOCADATA_H
#define ROOT_ACDPOCADATA_H

#include "TObject.h"
/** 
 * @class AcdPocaData
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

class AcdPocaData {
  
public:

  /// Default constructor.  
  AcdPocaData();
  
  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdPocaData(Float_t arcLength, Float_t doca, Float_t docaErr, 
              const TVector3& poca, const TVector3& pocaVector);

  /// Copy constructor
  AcdPocaData(const AcdPocaData& other);

  virtual ~AcdPocaData() {;}
  
  /// Assignment operator
  AcdPocaData& operator=(const AcdPocaData& other);
    
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
  inline Float_t getDocaErr() const {
    return m_docaErr;
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
  void set(Float_t arcLength, Float_t doca, Float_t docaErr, 
           const TVector3& poca, const TVector3& pocaVector);

  /// copy in the object
  void copy(const AcdPocaData& other);

  /// set individaul values
  void setArcLength(Float_t val) { m_arcLength = val; }
  void setDoca(Float_t val) { m_doca = val; }
  void setDocaErr(Float_t val) { m_docaErr = val; }
  void setPoca(const TVector3& val) { m_poca = val; }
  void setPocaVector(const TVector3& val) { m_pocaVector = val; }

  /// direct acces to members (to set them quickly)
  inline TVector3& accessPoca() { return m_poca; }
  inline TVector3& accessPocaVector() { return m_pocaVector; }

  // for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdPocaData&, const std::string & name = "" ) const ; // for tests

  /// reset all the values to their default (lowercase to avoid conflict with TObject::Clear())
  void clear();
  
  /// Print out this structure (lowercase to avoid conflict with TObject::Clear())
  void print(Option_t *option) const;
  
private:
  
  /// The arclength at which the poca occurs
  Float_t m_arcLength;

  /// The Distance of Closest Approach between track and tile or ribbon
  Float_t m_doca;                 
  
  /// The Error on the DOCA, this is the track error ellipsoid project along the doca direction
  Float_t m_docaErr;
  
  /// The Point of Closest Approach
  TVector3 m_poca;

  /// The vector from the POCA to the Detector Element edge
  TVector3 m_pocaVector;
  
  ClassDef(AcdPocaData,1)
  
};
#endif
