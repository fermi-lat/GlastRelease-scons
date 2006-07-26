
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

    void initXyz(Double_t x, Double_t y, Double_t z);

    UInt_t getSize() const { return m_taggerHitCol->GetEntriesFast();  }

    Double_t getPosition() const { return m_baricenterPosition; }

    Double_t getPulseHeight() const { return m_totalPulseHeight; };

    UInt_t getModuleId() const { return m_highestHit.getModuleId(); }

    UInt_t getLayerId() const { return m_highestHit.getLayerId(); }
    
    const commonRootData::TaggerHit& getHighestHit() const { return m_highestHit; };
    void setHighestHit(const commonRootData::TaggerHit& hit) { m_highestHit = hit; };
    Double_t getEta() const { return m_eta; }

    Double_t getX() const { return m_x; }
    Double_t getY() const { return m_y; }
    Double_t getZ() const { return m_z; }

    const TClonesArray* getTaggerHitCol() const { return m_taggerHitCol; }
    const commonRootData::TaggerHit* getTaggerHit(UInt_t ind) const;
    commonRootData::TaggerHit* addTaggerHit(UInt_t moduleId,UInt_t layerId,
                        UInt_t stripId, Double_t pulseHgt, Double_t sigma,
                        Bool_t isPedSub);


    void Fake( Int_t ievent, Float_t randNum );
    Bool_t CompareInRange( const TaggerCluster &, const std::string & name = "" )const;

private:

    Double_t m_baricenterPosition;
    Double_t m_totalPulseHeight;
    Double_t m_eta;
    Double32_t m_x;
    Double32_t m_y;
    Double32_t m_z;
    Int_t m_numTaggerHit;
    commonRootData::TaggerHit m_highestHit;

    TClonesArray *m_taggerHitCol;

    ClassDef(TaggerCluster,1) // ReconAncillary data beamtest 2006
};

} // end namespace
#endif
