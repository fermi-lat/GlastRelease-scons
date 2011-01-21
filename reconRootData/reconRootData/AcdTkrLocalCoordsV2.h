#ifndef ROOT_ACDTKRLOCALCOORDSV2_H
#define ROOT_ACDTKRLOCALCOORDSV2_H

#include "Rtypes.h"
#include <string>

#include <TMatrixDSym.h>
#include <TVector3.h>

/**
*  @class AcdTkrLocalCoords
*
*  @brief ROOT object which stores information about where a track crossed the plane of an ACD element.
*
*  This information is given in 2D.  The active distances are reported in the ACD element frame.
*  Both AcdTkrHitPoca and AcdTkrGapPoca inherit from this class.
*
*  The main access functions are:
*    - float getActiveX() , float  getActiveY()
*      - which return the active distance int the local frame of the ACD element
*        Positive values mean the track went into the volume in question, 
*        Negative values mean it missed the volume.
*        This magnitude is the distance to the edge of the volume.
*    - float getLocalXXCov() , float getLocalXYCov() , float getLocalYYCov()
*      - which return the terms of the covarience matrix projected into the place of the ACD element.
*    - float getCosTheta()
*      - which returns the angle of the track relative to detector plane
*    - float getPathLength()
*      - which returns the pathlength of track in detector element
*    - int getRegion()
*      - which returns a code telling which region of the tile was hit
*
*  
*  \author Eric Charles
*
* $Header$
*/

class AcdTkrLocalCoordsV2  {
  
public:

  AcdTkrLocalCoordsV2();

  AcdTkrLocalCoordsV2(Int_t volume, Float_t arcLength, Float_t cosTheta, 
                      const TVector3& global, 
                      const Float_t local[2], const Float_t active[2],
                      const TMatrixDSym& localCovProj, const TMatrixDSym& localCovProp);
  
  AcdTkrLocalCoordsV2(const AcdTkrLocalCoordsV2& other);
  
  virtual ~AcdTkrLocalCoordsV2() {};
  
  AcdTkrLocalCoordsV2& operator=(const AcdTkrLocalCoordsV2& other);
  
  /// Direct access to parameters
  
  /// Which volume
  inline Int_t  getVolume()             const {return m_volume; }
  
  /// Arclength of the track to the intersection
  inline float  getArclengthToPlane()    const {return m_arcLengthPlane; }

  /// Angle of track relative to detector plane
  inline float  getCosTheta()            const {return m_cosTheta; }

  // Location of hit in global coordinates
  inline const TVector3& getGlobalPosition() const {return m_global; };
     
  /// Location (and errors) of hit in tile coordinates
  inline float  getLocalX()             const {return m_local[0]; }
  inline float  getLocalY()             const {return m_local[1]; }

  /// Location (and errors) of hit in tile coordinates
  inline float  getActiveX()             const {return m_active[0]; }
  inline float  getActiveY()             const {return m_active[1]; }

  inline const TMatrixDSym& getLocalCovProj() const {return m_localCovProj; };
  inline const TMatrixDSym& getLocalCovProp() const {return m_localCovProp; };

  /// set everything at once
  void set(Int_t volume, Float_t arcLength, Float_t cosTheta, 
           const TVector3& global, 
           const Float_t local[2], const Float_t active[2],
           const TMatrixDSym& localCovProj, const TMatrixDSym& localCovProp);
  
  void copy(const AcdTkrLocalCoordsV2& other);
    
  // for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  
  Bool_t CompareInRange( const AcdTkrLocalCoordsV2 &, const std::string & name = "" ) const ; // for tests
  
  /// reset all the values to their default (lowercase to avoid conflict with TObject::Clear())
  void clear();
  
  /// Print out this structure (lowercase to avoid conflict with TObject::Clear())
  void print(Option_t *option="") const;
  
private:
  
  ///  Code that tells which part of tile was hit
  Int_t      m_volume;

  ///  Arclength to plane
  Float_t    m_arcLengthPlane;
 
  ///  Angle of track w.r.t. detector plane
  Float_t    m_cosTheta;

  ///  Position of the expected hit in Global Coords
  TVector3   m_global;

  ///  Position of the expected hit in Local Coords
  Float_t    m_local[2];  

  ///  Active distance of the hit
  Float_t    m_active[2];

  ///  Covariance terms in expected intersection
  TMatrixDSym m_localCovProj;
  
  ///  Covariance terms in expected intersection
  TMatrixDSym m_localCovProp;


  ClassDef(AcdTkrLocalCoordsV2,1)
    
};

#endif
