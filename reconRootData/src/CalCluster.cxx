#include "reconRootData/CalCluster.h"
#include <iostream>

ClassImp(CalCluster)

CalCluster::CalCluster() {
    Clear();
}

CalCluster::CalCluster(Double_t e, TVector3 p)
{ 
    Clear();
    m_energySum = e;
    m_energyCorrected = m_energySum;
    m_position = p;
}


void CalCluster::initialize(Double_t eLeak, Double_t eCorrect, Double_t rmsLong, Double_t rmsTrans, 
        const TVector3 &caldir, Double_t calTransvOffset) {
    
    m_leakEnergy = eLeak;
    m_energyCorrected = eCorrect;
    m_rmsLong = rmsLong;
    m_rmsTrans= rmsTrans;
    m_transvOffset = calTransvOffset;  
    m_direction = caldir;
}

void CalCluster::initialize(Double_t eLeak, Double_t eCorrect, std::vector<Double_t> eneLayer, 
        std::vector<TVector3> pLayer, std::vector<TVector3> rmsLayer, 
        Double_t rmsLong, Double_t rmsTrans, 
        const TVector3 &caldir, Double_t calTransvOffset) {
    
    m_leakEnergy = eLeak;
    m_energyCorrected = eCorrect;
    m_eneLayer = eneLayer;
    m_pLayer = pLayer;
    m_rmsLayer = rmsLayer;	
    m_rmsLong = rmsLong;
    m_rmsTrans= rmsTrans;
    m_transvOffset = calTransvOffset;  
    m_direction = caldir;
}

void CalCluster::initProfile(Double_t fitEnergy, Double_t ki2, Double_t fitStart,
                             Double_t fitAlpha, Double_t fitLambda) {
    m_fitEnergy = fitEnergy;
    m_profChisq = ki2;
    m_csiAlpha = fitAlpha;
    m_csiLambda = fitLambda;
    m_start  = fitStart;
    
}


void CalCluster::Clear(Option_t* /* option */) {
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

    m_eneLayer.clear();
    m_pLayer.clear();
    m_rmsLayer.clear();	


}

void CalCluster::Print(Option_t* /* option */) const {
    using namespace std;
    cout << "Energy: " << m_energySum << " Corrected: " << m_energyCorrected << endl;
    cout << "Pos: (" << m_position.X() << "," << m_position.Y() << "," 
        << m_position.Z() << ")" << endl;
    cout << "Dir: (" << m_direction.X() << "," << m_direction.Y() << "," 
        << m_direction.Z() << ")" << endl;
    cout << "Leak Energy: " << m_leakEnergy << "  RMS Long: " << m_rmsLong 
        << "  RMS Trans: " << m_rmsTrans << "  TransOffset: " << m_transvOffset
        << std::endl;
    cout << "Size of Energy Layer collection: " << m_eneLayer.size() << std::endl;
    cout << "Size of POS Layer collection: " << m_pLayer.size() << std::endl;
    cout << "Size of RMS Layer collection: " << m_rmsLayer.size() << std::endl;
}

