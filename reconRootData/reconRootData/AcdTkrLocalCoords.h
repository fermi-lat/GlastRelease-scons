#ifndef ROOT_ACDTKRLOCALCOORDS_H
#define ROOT_ACDTKRLOCALCOORDS_H

#include "Rtypes.h"
#include <string>

class TMatrixD;

/**
*  @class AcdTkrLocalCoords
*
*
*  @brief This class stores the data about the intersection of a track with an Acd Element in the local frame of that element
*  
*  \author Eric Charles
*
* $Header$
*/

class AcdTkrLocalCoords  {
  
public:

  AcdTkrLocalCoords();

  AcdTkrLocalCoords(const float position[2], float pathLength, float cosTheta, 
		    int region, const TMatrixD& covMatrix);
  
  AcdTkrLocalCoords(const AcdTkrLocalCoords& other);
  
  virtual ~AcdTkrLocalCoords() {};
  
  AcdTkrLocalCoords& operator=(const AcdTkrLocalCoords& other);
  
  /// Direct access to parameters
  
  /// Location (and errors) of hit in tile coordinates
  inline float  getActiveX()             const {return m_activeX; }
  inline float  getActiveY()             const {return m_activeY; }
  inline float  getLocalXXCov()          const {return m_XXCov; }
  inline float  getLocalYYCov()          const {return m_YYCov; }
  inline float  getLocalXYCov()          const {return m_XYCov; }
  
 /// Pathlength of track in detector element
  inline float  getPathLength()          const {return m_pathLength; }

  /// Angle of track relative to detector plane
  inline float  getCosTheta()            const {return m_cosTheta; }
     
  /// A code which tells which region of the tile was hit
  inline int    getRegion()              const {return m_region; }
  
  /// set everything at once
  void set(const float position[2], float pathLength, float cosTheta, 
	   int region, const TMatrixD& CovMatrix);
  
  void copy(const AcdTkrLocalCoords& other);
  
  // set the individual values (uncomment as needed)
  //inline void setPosition(const float position[2]) {
  //  m_X = position[0];
  //  m_Y = position[1];
  //};
  //inline void setRegion(int val) {m_region = val;}
  //inline void setCosTheta(float val) { m_cosTheta = val; }
  //inline void setPathLength(float val) { m_pathLength = val; }
  inline void setCovTerms(float XX, float YY, float XY) {
    m_XXCov = XX;
    m_YYCov = YY;
    m_XYCov = XY;      
  }
  
  // for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  
  Bool_t CompareInRange( const AcdTkrLocalCoords &, const std::string & name = "" ) const ; // for tests
  
  /// reset all the values to their default (lowercase to avoid conflict with TObject::Clear())
  void clear();
  
  /// Print out this structure (lowercase to avoid conflict with TObject::Clear())
  void print(Option_t *option="") const;
  
private:
  
  ///  X position of the expected hit in Local Coords
  Float_t    m_activeX;
  ///  Y position of the expected hit in Local Coords
  Float_t    m_activeY;
  
  ///  Pathlength of track through tile
  Float_t    m_pathLength;
  
  ///  Angle of track w.r.t. detector plane
  Float_t    m_cosTheta;

  ///  Covariance terms in expected LocalCoords
  Float_t    m_XXCov;        //  X Error squared  (x for Top, x for +-y planes, y for +- x plane) 
  Float_t    m_YYCov;        //  Y Error squared  (y for Top, z for +-x and +-y planes)
  Float_t    m_XYCov;        // correlation term of  X-Y error   
 
  ///  Code that tells which part of tile was hit
  Char_t     m_region;      
 
  ClassDef(AcdTkrLocalCoords,1)
    
};

#endif
