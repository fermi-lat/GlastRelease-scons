#ifndef ROOT_ACDTKRGAPPOCAV2_H
#define ROOT_ACDTKRGAPPOCAV2_H

#include "TObject.h"
#include "AcdTkrLocalCoordsV2.h"
#include "AcdPocaDataV2.h"

/** 
 * @class AcdTkrGapPocaV2
 * @brief ROOT object which stores information about the Point of Closest Approach (POCA) between 
 * an extrapolated track and a Gap in the ACD
 *  
 * Most of the structure of the object comes from the base classes AcdTkrLocalCoords and AcdPocaData
 * 
 * The class adds only enough information to define the involved in the POCA.
 *    - const AcdGapId& getId()  
 *      - which returns the ID of the gap
 *    - const int getTrackIndex()  
 *      - which returns the index of the track which did the hitting
 *
 * 
 * @author Eric Charles
 *
  */

#include "commonRootData/idents/AcdGapId.h"

class AcdTkrGapPocaV2 : public TObject, public AcdTkrLocalCoordsV2, public AcdPocaDataV2 {

public:

  /// Default constructor.  
  AcdTkrGapPocaV2();
  
  /// Copy constructor
  AcdTkrGapPocaV2(const AcdTkrGapPocaV2& params);

  /// Assignment operator
  AcdTkrGapPocaV2& operator=(const AcdTkrGapPocaV2& params);

  /// Return the AcdGapId of the Gap
  inline const AcdGapId& getGapId() const { return m_id; }
  
  /// Return the index of the associated track
  inline Int_t getTrackIndex() const {
    return m_trackIndex;
  }

  /// An estimator of the chance that a track in this gap would produce no signal
  inline Float_t vetoSigmaHit() const {
    return m_vetoSigmaHit;
  }
  
  /// An estimator of the chance that a track went into this gap
  inline Float_t vetoSigmaProj() const {
    return m_vetoSigmaProj;
  }
  
  /// An estimator of the chance that a track went into this gap
  inline Float_t vetoSigmaProp() const {
    return m_vetoSigmaProp;
  } 

  /// set all the values
  void initialize(const AcdGapId& acdId, int trackIndex,
                  Float_t vetoHit, Float_t vetoProj, Float_t vetoProp);
  
  /// reset all the values to their default
  virtual void Clear(Option_t*);
  
  /// Print out this structure
  virtual void Print(Option_t *option) const;
  
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdTkrGapPocaV2&, const std::string & name = "" ) const ; // for tests   

private:
  
  /// The ID of the Gap
  AcdGapId m_id;
  
  /// The index of the associated track
  Int_t m_trackIndex;
  
  /// An estimator of the chance that a track in this gap would produce no signal
  Float_t m_vetoSigmaHit;
  
  /// An estimator of the chance that a track went into this gap
  Float_t m_vetoSigmaProj;
  
    /// An estimator of the chance that a track went into this gap
  Float_t m_vetoSigmaProp;

  ClassDef(AcdTkrGapPocaV2,1)
  
};
#endif
