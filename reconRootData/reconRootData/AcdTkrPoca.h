#ifndef ROOT_ACDTKRPOCA_H
#define ROOT_ACDTKRPOCA_H

#include "TObject.h"
/** 
 * @class AcdTkrPoca
 * @brief Root object information about the Point of Closest Approach (POCA) between an extrapolated track
 *  and a hit Acd element (tile or ribbon).  This POCA is calculated in 3D.  The doca is defined to be positive 
 *  if the track goes inside the active distance and negative otherwise
 *  
 * This class should be a duplicate of Event::AcdTkrPoca
 * 
 * @author Eric Charles
 *
 * $Header$
 */

#include "commonRootData/idents/AcdId.h"
#include "TkrTrackParams.h"
#include "TVector3.h"

class AcdTkrPoca : public TObject {

public:

  /// what region of a tile does the poca occur in
  typedef enum {
    NONE_TILE = 0,
    HIT_TILE = 1,
    MISS_TILE = 2,      
    HIT_TOP, HIT_RIGHT, HIT_BOTTOM, HIT_LEFT,
    MISS_TOP, MISS_TOP_RIGHT, MISS_RIGHT, MISS_BOTTOM_RIGHT,
    MISS_BOTTOM, MISS_BOTTOM_LEFT, MISS_LEFT, MISS_TOP_LEFT } DocaRegionTile;
  
  /// what region of a ribbon does the poca occur in
  typedef enum {
    NONE_RIBBON = 0,
    HIT_RIBBON = 1,
    MISS_RIBBON = 2,      
    MISS_POSITIVE, MISS_NEGATIVE } DocaRegionRibbon;  
  
public:

  /// Default constructor.  
  AcdTkrPoca();
  
  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdTkrPoca(const AcdId& acdId, Int_t tkIndex, 
	     Float_t m_doca, Float_t m_docaErr, UInt_t docaRegion,
	     const TVector3& poca, const TkrTrackParams& paramsAtPoca);

  /// Copy constructor
  AcdTkrPoca(const AcdTkrPoca& params);

  /// Assignment operator
  AcdTkrPoca& operator=(const AcdTkrPoca& params);

  /// Return the AcdId of the hit tile or ribbon
  inline const AcdId& getId() const { return m_id; }
  
  /// Return the index of the associated track
  inline Int_t getTkrIndex() const {
    return m_tkIndex;
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
  
  /// Return a code showing which region of the tile or ribbon the poca occured in. 
  inline UInt_t getDocaRegion() const {
    return m_docaRegion;
  }
  
  /// Return the POCA (in global coordinates)
  const TVector3& getPoca() const {
    return m_poca;
  }
  
  /// Return the kalman propagated track parameters at the POCA
  const TkrTrackParams& getParamsAtPoca() const {
    return m_paramsAtPoca;
  }

  /// set all the values
  void initialize(const AcdId& acdId, Int_t tkIndex, 
		  Float_t doca, Float_t docaErr, UInt_t docaRegion,
		  const TVector3& poca, const TkrTrackParams& paramsAtPoca);

  /// reset all the values to their default
  virtual void Clear(Option_t*);
  
  /// Print out this structure
  virtual void Print(Option_t *option) const;
  
private:
  
  /// The ID of the hit tile
  AcdId m_id;
  
  /// The index of the associated track
  Int_t m_tkIndex;

  /// The Distance of Closest Approach between track and tile or ribbon
  Float_t m_doca;                 
  
  /// The Error on the DOCA, this is the track error ellipsoid project along the doca direction
  Float_t m_docaErr;
  
  /// The region relative to the tile in which the doca occured, see the enum above
  UInt_t m_docaRegion;
  
  /// The Point of Closest Approach
  TVector3 m_poca;
  
  /// The full track parameters at the POCA
  TkrTrackParams m_paramsAtPoca;

  ClassDef(AcdTkrPoca,1)
  
};
#endif
