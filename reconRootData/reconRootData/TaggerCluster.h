
#ifndef ROOT_TaggerCluster_H
#define ROOT_TaggerCluster_H

#include "TObject.h"
#include "TClonesArray.h"

#include "commonRootData/adf/TaggerHit.h"

/** @class TaggerCluster
 * @brief TaggerCluster data for beamtest 2006  
 * 
 * $Header$
*/
namespace reconRootData {

class TaggerCluster: public TObject {

public:
    
    TaggerCluster();
    virtual ~TaggerCluster();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void initialize(Double_t pos, Double_t pulseHeight, Double_t eta);

    UInt_t getSize() const { return m_taggerHitCol->GetEntriesFast();  }
    Double_t getPosition() const { return m_baricenterPosition; }
    Double_t getPulseHeight() const { return m_totalPulseHeight; };
    const commonRootData::TaggerHit& getHighestHit() const { return m_highestHit; };
    void setHighestHit(const commonRootData::TaggerHit& hit) { m_highestHit = hit; };
    Double_t getEta() const { return m_eta; }

    const TClonesArray* getTaggerHitCol() const { return m_taggerHitCol; }
    const commonRootData::TaggerHit* getTaggerHit(UInt_t ind) const;
    commonRootData::TaggerHit* addTaggerHit(UInt_t moduleId,UInt_t layerId,
                        UInt_t stripId,  UInt_t pulseHgt, Bool_t isPedSub);


    void Fake( Int_t ievent, Float_t randNum );
    Bool_t CompareInRange( const TaggerCluster &, const std::string & name = "" )const;

private:

    Double_t m_baricenterPosition;
    Double_t m_totalPulseHeight;
    Double_t m_eta;
    commonRootData::TaggerHit m_highestHit;
    Int_t m_numTaggerHit;

    TClonesArray *m_taggerHitCol;

    ClassDef(TaggerCluster,1) // ReconAncillary data beamtest 2006
};

} // end namespace
#endif
