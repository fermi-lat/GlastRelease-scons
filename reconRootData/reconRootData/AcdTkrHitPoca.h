#ifndef ROOT_ACDTKRHITPOCA_H
#define ROOT_ACDTKRHITPOCA_H

#include "TObject.h"
#include "AcdTkrLocalCoords.h"
#include "AcdPocaData.h"

/** 
 * @class AcdTkrHitPoca
 * @brief Root object information about the Point of Closest Approach (POCA) between an extrapolated track
 *  and a hit Acd element (tile or ribbon).  This POCA is calculated in 3D.  The doca is defined to be positive 
 *  if the track goes inside the active distance and negative otherwise
 *  
 * This class should be a duplicate of Event::AcdTkrHitPoca
 * 
 * @author Eric Charles
 *
 * $Header$
 */

#include "commonRootData/idents/AcdId.h"
#include "TkrTrackParams.h"
#include "TVector3.h"

class AcdTkrHitPoca : public TObject, public AcdTkrLocalCoords, public AcdPocaData {

public:

  /// Default constructor.  
  AcdTkrHitPoca();
  
  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdTkrHitPoca(const AcdId& acdId, int trackIndex,
		const AcdTkrLocalCoords& local,
		const AcdPocaData& pocaData);

  /// Copy constructor
  AcdTkrHitPoca(const AcdTkrHitPoca& params);

  /// Assignment operator
  AcdTkrHitPoca& operator=(const AcdTkrHitPoca& params);

  /// Return the AcdId of the hit tile or ribbon
  inline const AcdId& getId() const { return m_id; }
  
  /// Return the index of the associated track
  inline Int_t getTrackIndex() const {
    return m_trackIndex;
  }

  /// set all the values
  void set(const AcdId& acdId, int trackIndex,
	   const AcdTkrLocalCoords& local,
	   const AcdPocaData& pocaData);
  
  /// reset all the values to their default
  virtual void Clear(Option_t*);
  
  /// Print out this structure
  virtual void Print(Option_t *option) const;
  
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdTkrHitPoca&, const std::string & name = "" ) const ; // for tests   

private:
  
  /// The ID of the hit tile
  AcdId m_id;
  
  /// The index of the associated track
  Int_t m_trackIndex;

  ClassDef(AcdTkrHitPoca,1)
  
};
#endif
