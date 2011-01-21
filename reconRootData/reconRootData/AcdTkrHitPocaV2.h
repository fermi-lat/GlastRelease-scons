#ifndef ROOT_ACDTKRHITPOCAV2_H
#define ROOT_ACDTKRHITPOCAV2_H

#include "TObject.h"
#include "AcdTkrLocalCoordsV2.h"
#include "AcdPocaDataV2.h"

/** 
 * @class AcdTkrHitPocaV2
 * @brief ROOT object which stores information about the Point of Closest Approach (POCA) between an 
 * extrapolated track and a hit Acd element (tile or ribbon).  
 *  
 * Most of the structure of the object comes from the base classes AcdTkrLocalCoords and AcdPocaData
 * 
 * The class adds only enough information to define the involved in the POCA.
 *    - const idents::AcdId& getId()  
 *      - which returns the ID of the hit element
 *    - const int getTrackIndex()  
 *      - which returns the index of the track which did the hitting
 *
 * @author Eric Charles
 *
 * $Header$
 */

#include "commonRootData/idents/AcdId.h"
#include "TkrTrackParams.h"
#include "TVector3.h"

class AcdTkrHitPocaV2 : public TObject, public AcdTkrLocalCoordsV2, public AcdPocaDataV2 {

public:

  /// Default constructor.  
  AcdTkrHitPocaV2();
  
  /// Copy constructor
  AcdTkrHitPocaV2(const AcdTkrHitPocaV2& params);

  /// Assignment operator
  AcdTkrHitPocaV2& operator=(const AcdTkrHitPocaV2& params);

  /// Return the AcdId of the hit tile or ribbon
  inline const AcdId& getId() const { return m_id; }
  
  /// Return the index of the associated track
  inline Int_t getTrackIndex() const {
    return m_trackIndex; 
  }  
  
  /// Return the mips associated with PMT A
  inline Float_t mipsPmtA() const { 
    return m_mips[0];
  }
  
  /// Return the mips associated with PMT B
  inline Float_t mipsPmtB() const { 
    return m_mips[1];
  }
 
  /// An estimator of the number of sigma needed for this track to be a true MIP signal
  inline Float_t vetoSigmaHit() const {
    return m_vetoSigmaHit;
  }
  
  /// An estimator of the number of sigma needed for this track to hit this element
  inline Float_t vetoSigmaProj() const {
    return m_vetoSigmaProj;
  }
  
  /// An estimator of the number of sigma needed for this track to hit this element
  inline Float_t vetoSigmaProp() const {
    return m_vetoSigmaProp;
  } 

  /// set all the values
  void set(const AcdId& acdId, int trackIndex,
           const float mips[2],
           Float_t vetoHit, Float_t vetoProj, Float_t vetoProp);
  
  /// reset all the values to their default
  virtual void Clear(Option_t*);
  
  /// Print out this structure
  virtual void Print(Option_t *option) const;
  
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdTkrHitPocaV2&, const std::string & name = "" ) const ; // for tests   

private:
  
  /// The ID of the hit tile
  AcdId m_id;
  
  /// The index of the associated track
  Int_t m_trackIndex;
         
  /// The mip values associated with the A pmt
  Float_t m_mips[2];

  ///  An estimator of the number of sigma needed for this hit to be a true MIP signal
  Float_t m_vetoSigmaHit;
  
  ///  An estimator of the number of sigma needed for this track to hit this element
  Float_t m_vetoSigmaProj;
  
  ///  An estimator of the number of sigma needed for this track to hit this element
  Float_t m_vetoSigmaProp;

  ClassDef(AcdTkrHitPocaV2,1)
  
};
#endif
