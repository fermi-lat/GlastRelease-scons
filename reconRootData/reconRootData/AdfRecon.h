
#ifndef ROOT_AdfRecon_H
#define ROOT_AdfRecon_H

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"

#include "TaggerCluster.h"
#include "commonRootData/adf/QdcHit.h"
#include "commonRootData/adf/ScalerHit.h"

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

    void initEventNumber(UInt_t evtNum) { m_eventNumber = evtNum;}
    void initSpillNumber(UInt_t spillNum) { m_spillNumber = spillNum; }
    void initNmodule(UInt_t n) { m_nModule = n; }

    UInt_t getEventNumber() const { return m_eventNumber; }
    UInt_t getSpillNumber() const { return m_spillNumber; }
    UInt_t getNmodule() const { return m_nModule; }

    void initNumHighestCluster(UInt_t n) {m_numHighestCluster = n;}
    UInt_t getNumHighestCluster() const { return m_numHighestCluster; }

    void initXyz(const Double_t *x, const Double_t *y, const Double_t *z, UInt_t n);

    Double_t getX(Int_t module) const { 
        if ((module >= 0) && (module < m_nModule))
            return m_x[module]; 
        else return 0.0; }
    Double_t getY(Int_t module) const { 
        if((module >= 0) && (module < m_nModule))
           return m_y[module];
        else return 0.0;  }
    Double_t getZ(Int_t module) const { 
        if((module >= 0) && (module < m_nModule))
            return m_z[module]; 
        else return 0.0; }

    void init(Double_t px, Double_t py, Double_t pz) {
        m_px = px; m_py = py; m_pz = pz; }

    Double_t getPx() const { return m_px; }
    Double_t getPy() const { return m_py; }
    Double_t getPz() const { return m_pz; }

    void initPhi(Double_t in, Double_t out) {
        m_phiIn = in; m_phiOut = out; }

    Double_t getPhiIn() const { return m_phiIn; }
    Double_t getPhiOut() const { return m_phiOut; }

    void initThetaPhi(Double_t t, Double_t p) {
        m_theta = t; m_dphi = p; }

    Double_t getTheta() const { return m_theta; }
    Double_t getDeltaPhi() const { return m_dphi; }

    void initEnergy(Double_t eRec, Double_t eCorr) {
        m_eRec = eRec; m_eCorr = eCorr; }
    Double_t getReconstructedEnergy() const { return m_eRec; }
    Double_t getCorrectedEnergy() const { return m_eCorr; }

    const TObjArray* getTaggerClusterCol() const { return m_taggerClusterCol; }
    const TaggerCluster* getTaggerCluster(UInt_t ind) const;
    void addTaggerCluster(TaggerCluster *clus);


    const TClonesArray* getQdcHitCol() const { return m_qdcHitCol; }
    const commonRootData::QdcHit* getQdcHit(UInt_t ind) const;
    commonRootData::QdcHit* addQdcHit(UInt_t channel, Double_t pulseHgt, UInt_t module, Double_t sigma, Bool_t isPedSub);

    const TClonesArray* getScalerHitCol() const { return m_scalerHitCol; }
    const commonRootData::ScalerHit* getScalerHit(UInt_t ind) const;
    commonRootData::ScalerHit* addScalerHit(UInt_t chan, UInt_t val);

    void Fake( Int_t ievent, Float_t randNum );
    Bool_t CompareInRange( const AdfRecon &, const std::string & name = "" )const;

private:

    UInt_t m_eventNumber;
    UInt_t m_spillNumber;

    UInt_t m_numHighestCluster;
    Int_t m_numQdcHit;
    Int_t m_numScalerHit;

    Int_t m_nModule;
    Double32_t *m_x; //[m_nModule]
    Double32_t *m_y; //[m_nModule]
    Double32_t *m_z; //[m_nModule]

    Double32_t m_px;
    Double32_t m_py;
    Double32_t m_pz;
    Double32_t m_eRec;
    Double32_t m_eCorr;
    Double32_t m_phiIn;
    Double32_t m_phiOut;
    Double32_t m_theta;
    Double32_t m_dphi;

    TObjArray *m_taggerClusterCol;

    TClonesArray *m_qdcHitCol;
    TClonesArray *m_scalerHitCol;

    ClassDef(AdfRecon,1) // ReconAncillary data beamtest 2006
};

} // end namespace
#endif
