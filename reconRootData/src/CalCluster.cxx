#include "reconRootData/CalCluster.h"
#include <iostream>

ClassImp(CalCluster)

CalCluster::CalCluster() {
    Clear();
}

void CalCluster::initialize(Double_t eLeak, std::vector<Double_t> eneLayer, 
        std::vector<TVector3> pLayer, std::vector<TVector3> rmsLayer, 
        Double_t rmsLong, Double_t rmsTrans, 
        const TVector3 &caldir, Double_t calTransvOffset) {
    
    m_eneLayer = eneLayer;
    m_leakEnergy = eLeak;
    m_pLayer = pLayer;
    m_rmsLayer = rmsLayer;	
    m_rmsLong = rmsLong;
    m_rmsTrans= rmsTrans;
    m_transvOffset = calTransvOffset;    
}

void CalCluster::initProfile(Double_t fitEnergy, Double_t ki2, Double_t fitStart,
                             Double_t fitAlpha, Double_t fitLambda) {
    m_fitEnergy = fitEnergy;
    m_profChisq = ki2;
    m_csiAlpha = fitStart;
    m_csiLambda = fitAlpha;
    m_start  = fitLambda;
    
}


void CalCluster::Clear(Option_t *option) {
    m_energySum = 0.;
    m_leakEnergy = 0.;
    m_energyCorrected = 0.;
    m_rmsLong = 0.0;
    m_rmsTrans = 0.0;
    m_transvOffset = 0.0;
    m_fitEnergy = 0.0;
    m_profChisq = 0.0;
    m_csiAlpha = 0.0;
    m_csiLambda = 0.0;
    m_start = 0.0;
    
    m_position = TVector3(0.,0.,0.);
    m_direction = TVector3(0.,0.,0);

}

void CalCluster::Print(Option_t *option) const {
    using namespace std;
    cout << "Energy: " << m_energySum << " Corrected: " << m_energyCorrected << endl;
    cout << "Pos: (" << m_position.X() << "," << m_position.Y() << "," 
        << m_position.Z() << ")" << endl;
    cout << "Dir: (" << m_direction.X() << "," << m_direction.Y() << "," 
        << m_direction.Z() << ")" << endl;
}

