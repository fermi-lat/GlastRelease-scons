#ifndef CALCLUSTER_H
#define CALCLUSTER_H

// Root Includes
#include "TObject.h"
#include "TVector3.h"

//! Root object that stores a cluster of CsI detections
/*! Provides access to the positions and energy depositions for a
cluster of CsI logs */


#include <vector>
#ifdef WIN32
using namespace std;
#endif


class CalCluster : public TObject
{
    
private:
    
    //! Total measured energy in the calorimeter
    Double_t m_energySum;
    //! Leakage corrected energy using correlation method ( for E> several GeV)
    Double_t m_leakEnergy;
    //! corrected energy not used ( yet )
    Double_t m_energyCorrected;
    //! Energy per layer in MeV
    std::vector<Double_t> m_eneLayer;
    //! Barycenter position in each layer
    std::vector<TVector3> m_pLayer;
    //! RMS of energy deposition in each layer
    std::vector<TVector3> m_rmsLayer;
    //! RMS of longitudinal position measurement
    Double_t m_rmslong;
    //! RMS of transverse position measurement
    Double_t m_rmstrans;
    //! Transvers offset of calorimeter position measurement
    Double_t m_transvOffset;
    
    //! fitted energy ( for E>10 GeV)
    Double_t m_fitEnergy;
    //! Chisquare of the fit ( not a real Chisquare)
    Double_t m_ProfChisq;
    //! Alpha parameter used in the fit
    Double_t m_CsiAlpha;
    //! Lambda parameter used in the fit 
    Double_t m_CsiLambda;
    //! Fitted starting point of the shower (physical meaning is not clear)
    Double_t m_start;
    
    TVector3 m_position;
    TVector3 m_direction;
    
public:
 
  
    //! constructor
    CalCluster() {};
    CalCluster(Double_t e, TVector3 p);
    
    //! Destructor
    ~CalCluster() {};
    
    void setDirection(TVector3 v)   {m_direction = v;}
    
    //! Set energy corrected
    /*! not used for the moment
    *  Energy sum is stored */
    void setEnergyCorrected(Double_t e) {m_energyCorrected = e;}
    //! Set energy per layer
    void setEneLayer(std::vector<Double_t> v){m_eneLayer = v;}
    //! Set barycenter position for each layer
    void setPosLayer(std::vector<TVector3> v){m_pLayer = v;}
    //! Set rms of energy deposition for each layer
    void setRmsLayer(std::vector<TVector3> v){m_rmsLayer = v;}
    //! Set Longitudinal RMS
    void setRmsLong(Double_t r) {m_rmslong=r;}
    //! Set transverse RMS
    void setRmsTrans(Double_t r) {m_rmstrans=r;}
    //! Set energy corrected via CalClustersAlg::Leak()
    void setEneLeak(Double_t e) {m_leakEnergy = e;}
    //! Set fitted energy form CalClustersAlg::Profile()
    void setFitEnergy(Double_t e) { m_fitEnergy = e;}
    //! Set chi square of profile fitting
    void setProfChisq(Double_t k) { m_ProfChisq = k;}
    //! Set alpha parameter used in the fit
    void setCsiAlpha(Double_t a) { m_CsiAlpha =a;}
    //! Set lambda parameter used in the fit
    void setCsiLambda(Double_t l) { m_CsiLambda = l;}
    //! Set the fitted starting point
    void setCsiStart(Double_t s) { m_start = s;}
    //! Set the transverse offset of calorimeter position measurement
    void setTransvOffset (Double_t offset) {m_transvOffset = offset;}
    
    // access
    Double_t energySum()        const {return m_energySum;}
    Double_t energyLeak()       const {return m_leakEnergy;}
    Double_t energyCorrected()  const {return m_energyCorrected;}
    Double_t getEneLayer(int i) const {return m_eneLayer[i];}
    const TVector3& getPosLayer(int i) const {return m_pLayer[i];}
    const std::vector<Double_t>& getEneLayer() const {return m_eneLayer;}
    const std::vector<TVector3>& getPosLayer() const {return m_pLayer;}
    const std::vector<TVector3>& getRmsLayer() const {return m_rmsLayer;}
    Double_t getRmsLong()		  const {return m_rmslong;}
    Double_t getRmsTrans()	  const {return m_rmstrans;}
    Double_t getTransvOffset()  const {return m_transvOffset;}
    
    TVector3 position()          const {return m_position;}
    TVector3 direction()        const {return m_direction;}
    Double_t getFitEnergy()	  const {return m_fitEnergy;}
    Double_t getProfChisq()	  const {return m_ProfChisq;}
    Double_t getCsiAlpha()	  const {return m_CsiAlpha;}
    Double_t getCsiLambda()	  const {return m_CsiLambda;}
    Double_t getCsiStart()	  const {return m_start;}
    // operations
    void writeOut() const;
    
protected:
    
    virtual void ini();
    
    ClassDef(CalCluster,3)
};


#endif 


