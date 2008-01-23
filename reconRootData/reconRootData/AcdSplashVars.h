#ifndef ROOT_ACDBEAMVARS_H
#define ROOT_ACDBEAMVARS_H

#include "TObject.h"

/** 
 * @class AcdSplashVars
 * @brief A ROOT object which stores information about the geometry of potential backsplash into the ACD
 *
 * The data in this class are calculated by extrapolating a track to the CAL entry point
 * and then looking at the solid angle occupied by a tile w.r.t. that point.
 *
 *
 *  The main access functions are:
 *    - const idents::AcdId& getId()  
 *      - which returns the ID of the hit element
 *    - const int getTrackIndex()  
 *      - which returns the index of the track which did the hitting
 *    - const TVector3& calEntryPoint() 
 *      - which returns the point the track enters the calorimeter
 *    - const TVector3& calEntryVector() 
 *      - which returns the track vector at the point the track enters the calorimeter
 *    - float tileSolidAngle
 *      - which returns the total solid angle of the tile, seen from the track entry point
 *    - float weightedTrackAngle
 *      - which returns the average of the angle between the track vector and the vector from the CAL to the tile
 *    - float weightedPathlength
 *      - which returns the average of the pathlength inside the tile along the path from the CAL to the tile
 * 
 * @author Eric Charles
 *
 * $Header$
 */

#include "commonRootData/idents/AcdId.h"
#include "TVector3.h"

class AcdSplashVars : public TObject {

public:

  /// Default constructor.  
  AcdSplashVars();
  
  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdSplashVars(const AcdId& acdId, int trackIndex, 
	      const TVector3& calEntryPoint, const TVector3& calEntryVector,
	      const Float_t& tileSolidAngle, const Float_t& weightedTrackAngle,
	      const Float_t& weightedPathlength);

  /// Copy constructor
  AcdSplashVars(const AcdSplashVars& params);

  /// Assignment operator
  AcdSplashVars& operator=(const AcdSplashVars& params);

  /// Return the AcdId of the hit tile or ribbon
  inline const AcdId& getId() const { return m_id; }
  
  /// Return the index of the associated track
  inline Int_t getTrackIndex() const {
    return m_trackIndex;
  }

  /// return the Vector from the point the track enters the calorimeter to the tile center
  inline const TVector3& calEntryPoint() const { return m_calEntryPoint; }
  
  /// return the track vector at the point the track enters the calorimeter
  inline const TVector3& calEntryVector() const { return m_calEntryVector; }
  
  /// return total solid angle of the tile, seen from the track entry point
  inline const Float_t&  tileSolidAngle() const { return m_tileSolidAngle; }

  /// return the average of the angle between the reconstructed track and the line         
  /// connecting the track entry point in the calorimeter at that point
  /// (weighted by the solid angle of the element)                         
  inline const Float_t&  weightedTrackAngle() const { return m_weightedTrackAngle; }

  /// return the average of the pathlength inside the tile along the path from the 
  /// track entry point in the calorimeter to the tile
   /// (weighted by the solid angle of the element)     
  inline const Float_t&  weightedPathlength() const { return m_weightedPathlength; }
  
  /// set all the values
  void set(const AcdId& acdId, int trackIndex, 
	   const TVector3& calEntryPoint, const TVector3& calEntryVector,
	   const Float_t& tileSolidAngle, const Float_t& weightedTrackAngle,
	   const Float_t& weightedPathlength);
  
  /// reset all the values to their default
  virtual void Clear(Option_t*);
  
  /// Print out this structure
  virtual void Print(Option_t *option) const;
  
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdSplashVars&, const std::string & name = "" ) const ; // for tests   

private:
  
  /// The ID of the hit tile
  AcdId m_id;
  
  /// The index of the associated track
  Int_t m_trackIndex;

  /// The Vector from the point the track enters the calorimeter to the tile center
  TVector3 m_calEntryPoint;

  /// The track vector at the point the track enters the calorimeter
  TVector3 m_calEntryVector;
  
  /// The total solid angle of the tile, seen from the track entry point
  Float_t m_tileSolidAngle;

  /// The average of the angle between the reconstructed track and the line         
  /// connecting the track entry point in the calorimeter at that point
  /// (weighted by the solid angle of the element)                         
  Float_t m_weightedTrackAngle;

  /// the average of the pathlength inside the tile along the path from the 
  /// track entry point in the calorimeter to the tile
   /// (weighted by the solid angle of the element)     
  Float_t m_weightedPathlength;
  
  ClassDef(AcdSplashVars,1)
  
};
#endif
