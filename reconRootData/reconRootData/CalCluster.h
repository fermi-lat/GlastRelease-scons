#ifndef ROOT_CALCLUSTER_H
#define ROOT_CALCLUSTER_H

#include "TObject.h"
#include "TVector3.h"

#include <vector>
#ifdef WIN32
using namespace std;
#endif

/** @class CalCluster
 * @brief Root object that stores a cluster of CsI detections
 *
 * Provides access to the positions and energy depositions for a
 * cluster of CsI logs 
 *
 * @author A. Chekhtman
 *
 * $Header$
 */
class CalCluster : public TObject
{ 
public:
    
    CalCluster();
    
    CalCluster(Double_t e, TVector3 p);

    virtual ~CalCluster() {};

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;
    
    void initialize(Double_t eLeak, std::vector<Double_t> eneLayer, 
        std::vector<TVector3> pLayer, std::vector<TVector3> rmsLayer, 
        Double_t rmsLong, Double_t rmsTrans, 
        const TVector3 &caldir, Double_t calTransvOffset);

    void initialize(Double_t eLeak, Double_t rmsLong, Double_t rmsTrans, 
        const TVector3 &caldir, Double_t calTransvOffset);

    void initProfile(Double_t fitEnergy, Double_t ki2, Double_t fitStart,
        Double_t fitAlpha, Double_t fitLambda);

    void addEnergy(Double_t e) { m_eneLayer.push_back(e); };
    void addPos(const TVector3& pos) { m_pLayer.push_back(pos); };
    void addRms(const TVector3& rms) { m_rmsLayer.push_back(rms); };
    
    Double_t getEnergySum() const {return m_energySum;}
    Double_t getEnergyLeak() const {return m_leakEnergy;}
    Double_t getEnergyCorrected() const {return m_energyCorrected;}
    Double_t getEneLayer(int i) const {return m_eneLayer[i];}
    UInt_t getNumEneLayer() const { return m_eneLayer.size();}
    const TVector3& getPosLayer(int i) const {return m_pLayer[i];}
    UInt_t getNumPosLayer() const { return m_pLayer.size(); }
    const TVector3& getRmsLayer(int i) const { return m_rmsLayer[i];}
    UInt_t getNumRmsLayer() const { return m_rmsLayer.size(); }
    const std::vector<Double_t>& getEneLayer() const {return m_eneLayer;}
    const std::vector<TVector3>& getPosLayer() const {return m_pLayer;}
    const std::vector<TVector3>& getRmsLayer() const {return m_rmsLayer;}
    Double_t getRmsLong() const {return m_rmsLong;}
    Double_t getRmsTrans() const {return m_rmsTrans;}
    Double_t getTransvOffset() const {return m_transvOffset;}
    
    const TVector3& getPosition() const {return m_position;}
    const TVector3& getDirection() const {return m_direction;}
    Double_t getFitEnergy()	  const {return m_fitEnergy;}
    Double_t getProfChisq()	  const {return m_profChisq;}
    Double_t getCsiAlpha()	  const {return m_csiAlpha;}
    Double_t getCsiLambda()	  const {return m_csiLambda;}
    Double_t getCsiStart()	  const {return m_start;}        
    
private:
    
    /// Total measured energy in the calorimeter
    Double_t m_energySum;
    /// Leakage corrected energy using correlation method ( for E> several GeV)
    Double_t m_leakEnergy;
    /// corrected energy not used ( yet )
    Double_t m_energyCorrected;
    /// RMS of longitudinal position measurement
    Double_t m_rmsLong;
    /// RMS of transverse position measurement
    Double_t m_rmsTrans;
    /// Transvers offset of calorimeter position measurement
    Double_t m_transvOffset;
    
    /// fitted energy ( for E>10 GeV)
    Double_t m_fitEnergy;
    /// Chisquare of the fit ( not a real Chisquare)
    Double_t m_profChisq;
    /// Alpha parameter used in the fit
    Double_t m_csiAlpha;
    /// Lambda parameter used in the fit 
    Double_t m_csiLambda;
    /// Fitted starting point of the shower (physical meaning is not clear)
    Double_t m_start;
    
    TVector3 m_position;
    TVector3 m_direction;
    
    /// Energy per layer in MeV
    std::vector<Double_t> m_eneLayer;
    /// Barycenter position in each layer
    std::vector<TVector3> m_pLayer;
    /// RMS of energy deposition in each layer
    std::vector<TVector3> m_rmsLayer;

    ClassDef(CalCluster,3)
};


#endif 


