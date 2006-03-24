#ifndef ROOT_ACDTKRINTERSECTION_H
#define ROOT_ACDTKRINTERSECTION_H

#include "TObject.h"

#include "commonRootData/idents/AcdId.h"
#include "TVector3.h"

#include "TMatrixD.h"

/** 
 * @class AcdTkrIntersection
 * @brief Root object containing information about the expected intersection of a track with the ACD
 * 
 * This class should be a duplicate of Event::AcdTkrIntersection
 * 
 * @author Eric Charles
 *
 * $Header$
 */

class AcdTkrIntersection : public TObject
{
public:

    AcdTkrIntersection();

    AcdTkrIntersection(const AcdId& tileId, int trackIndex, 
		       const TVector3& globalPosition, 
		       const Double_t localPosition[2], const TMatrixD& localCovMatrix,
		       Double_t arcLengthToIntersection, Double_t pathLengthInTile, UChar_t tileHit);

    AcdTkrIntersection(const AcdTkrIntersection& params);

    virtual ~AcdTkrIntersection() {};

    /// resets all member variables to default values
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const AcdTkrIntersection &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    AcdTkrIntersection & operator=(const AcdTkrIntersection& old);

    /// Provide an initialization method (in case of default constructor)
    void initialize(const AcdId& tileId, Int_t trackIndex, 
		    const TVector3& globalPostion, 
		    const Double_t localPosition[2], const TMatrixD& localCovMatrix,
		    Double_t arcLengthToIntersection, Double_t pathLengthInTile, UChar_t tileHit);
    void initialize(const AcdTkrIntersection& params);

    /// Direct access to parameters

    // Which tile should be hit
    inline const AcdId&    getTileId()              const {return m_tileId; };
    inline Int_t           getTrackIndex()          const {return m_trackIndex; };

    // Location of hit in global coordinates
    inline const TVector3& getGlobalPosition()      const {return m_location;    };
    inline TVector3&       getGlobalPosition()            {return m_location;    };
    
    // Location (and errors) of hit in tile coordinates
    inline Double_t  getLocalX()              const {return m_localX; };
    inline Double_t  getLocalY()              const {return m_localY; };
    inline Double_t  getLocalXXCov()          const {return m_localXXCov; };
    inline Double_t  getLocalYYCov()          const {return m_localYYCov; };
    inline Double_t  getLocalXYCov()          const {return m_localXYCov; };

    // Distance along track from first hit to tile intersection
    inline Double_t  getArcLengthToIntersection() const { return m_arcLengthToIntersection; } ;
    // Path length of track through tile
    inline Double_t  getPathLengthInTile()    const { return m_pathlengthInTile; } ;

    // is this tile hit
    inline UChar_t tileHit() const { return m_tileHit; };

private:
  
    /// ID of hit tile
    AcdId      m_tileId;

    /// Index of the related track
    Int_t         m_trackIndex;

    /// Global position of expected hit 
    TVector3      m_location;
    ///  X Position of the expected hit in Tile Coords
    Double32_t    m_localX;
    ///  Y Position of the expected hit in Tile Coords
    Double32_t    m_localY;
    
    ///  Covariance terms in expected intersection
    Double32_t    m_localXXCov;        // local X Error squared  (x for Top, x for +-y planes, y for +- x plane) 
    Double32_t    m_localYYCov;        // local Y Error squared  (y for Top, z for +-x and +-y planes)
    Double32_t    m_localXYCov;        // correlation term of local X-Y error   

    ///  Distance from first hit to intersection
    Double32_t    m_arcLengthToIntersection;
  
    ///  Pathlength through the ACD tile
    Double32_t    m_pathlengthInTile;

    ///  Mask to say if the tile has a signal, and if it fired discriminators
    UChar_t       m_tileHit;

    ClassDef(AcdTkrIntersection,1)

};





#endif
