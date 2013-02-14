#ifndef ROOT_ACDASSOC_H
#define ROOT_ACDASSOC_H

#include "TObject.h"

/** 
 * @class AcdAssoc
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
 * @author Alex Drlica-Wagner
 *
 * $Header$
 */

#include "commonRootData/idents/AcdId.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TClonesArray.h"
#include "AcdTkrPointV2.h"


class AcdTkrHitPocaV2;
class AcdTkrGapPocaV2;

class AcdAssoc : public TObject {

public:

  /// Default constructor.  
  AcdAssoc();
  
  /// Copy constructor
  AcdAssoc(const AcdAssoc& params);

  /// Assignment operator
  AcdAssoc& operator=(const AcdAssoc& params);

  /// Return the index of the associated track
  inline Int_t getTrackIndex() const {
    return m_trackIndex; 
  }
  
  inline Bool_t getUpward() const { return m_upward; }
  
  inline Float_t getEnergy() const { return m_energy; }
  
  inline const TVector3& getStart() const { return m_start; }
  
  inline const TVector3& getDir() const { return m_dir; }
  
  inline Float_t getArcLength() const { return m_arcLength; }
  
  inline int getTkrSSDVeto() const { return m_tkrSSDVeto; }
  
  inline Float_t getCornerDoca() const { return m_cornerDoca; }

  inline const TMatrixDSym& getCovStart() const { return m_cov_start; }
  
  inline const TMatrixDSym& getCovEnd() const { return m_cov_end; }

  UInt_t nAcdHitPoca() const { return m_hitPocae.GetEntriesFast(); }

  UInt_t nAcdGapPoca() const { return m_gapPocae.GetEntriesFast(); }

  const AcdTkrHitPocaV2* getHitPoca(UInt_t i = 0) const;

  const AcdTkrGapPocaV2* getGapPoca(UInt_t i = 0) const;

  const AcdTkrPointV2& getPoint() const { return m_point; }
 
  AcdTkrHitPocaV2* addHitPoca(const AcdTkrHitPocaV2& poca);

  AcdTkrGapPocaV2* addGapPoca(const AcdTkrGapPocaV2& poca);

  void setPoint(const AcdTkrPointV2& point);

  /// set all the values
  void set(Int_t index, Bool_t up, Float_t energy, 
           const TVector3& start, const TVector3& dir, Float_t arcLength,
           const TMatrixDSym& covStart, const TMatrixDSym& covEnd,
           Int_t tkrSSDVeto, Float_t cornerDoca);
  
  /// reset all the values to their default
  virtual void Clear(Option_t*);
  
  /// Print out this structure
  virtual void Print(Option_t *option) const;
  
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const AcdAssoc&, const std::string & name = "" ) const ; // for tests   

private:
  

  Int_t           m_trackIndex;
  
  Bool_t          m_upward;
  
  Float_t         m_energy;
  
  TVector3        m_start;

  TVector3        m_dir;
  
  Float_t         m_arcLength;
  
  TMatrixDSym     m_cov_start;
  
  TMatrixDSym     m_cov_end;     
  
  Int_t           m_tkrSSDVeto;

  Float_t         m_cornerDoca;

  TClonesArray    m_hitPocae;
  
  TClonesArray    m_gapPocae;

  AcdTkrPointV2   m_point;  

  ClassDef(AcdAssoc,1)
  
};
#endif
