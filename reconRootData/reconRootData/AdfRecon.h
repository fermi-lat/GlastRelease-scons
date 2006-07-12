
#ifndef ROOT_AdfRecon_H
#define ROOT_AdfRecon_H

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"

#include "TaggerCluster.h"
#include "commonRootData/adf/QdcHit.h"

/** @class AdfRecon
 * @brief The digitization ancillary data for beamtest 2006  
 * 
 * $Header$
*/

namespace reconRootData {

class AdfRecon: public TObject {

public:
    
    AdfRecon();
    virtual ~AdfRecon();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void setEventNumber(UInt_t evtNum) { m_eventNumber = evtNum;}
    void setSpillNumber(UInt_t spillNum) { m_spillNumber = spillNum; }

    UInt_t getEventNumber() const { return m_eventNumber; }
    UInt_t getSpillNumber() const { return m_spillNumber; }

    const TObjArray* getTaggerClusterCol() const { return m_taggerClusterCol; }
    const TaggerCluster* getTaggerCluster(UInt_t ind) const;
    void addTaggerCluster(TaggerCluster *clus);


    const TClonesArray* getQdcHitCol() const { return m_qdcHitCol; }
    const commonRootData::QdcHit* getQdcHit(UInt_t ind) const;
    commonRootData::QdcHit* addQdcHit(UInt_t channel, UInt_t pulseHgt, Bool_t isPedSub);


    void Fake( Int_t ievent, Float_t randNum );
    Bool_t CompareInRange( const AdfRecon &, const std::string & name = "" )const;

private:

    UInt_t m_eventNumber;
    UInt_t m_spillNumber;

    Int_t m_numQdcHit;

    TObjArray *m_taggerClusterCol;

    TClonesArray *m_qdcHitCol;

    ClassDef(AdfRecon,1) // ReconAncillary data beamtest 2006
};

} // end namespace
#endif
