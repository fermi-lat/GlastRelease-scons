#include "reconRootData/CalCluster.h"


ClassImp(CalCluster)


CalCluster::CalCluster (Int_t ind)
:m_xposLayer(8),m_yposLayer(8),m_eneLayer(8)
{
    m_position = new TVector3(0,0,0);
    m_direction = new TVector3(0,0,0);
    m_id = ind;
    m_Esum = 0.0f;
    m_Ecorr = 0.0f;
	m_CsiAlpha=0.0f;
	m_CsiLambda=0.0f;
	m_start=0.0f;
	m_ProfChisq=0.0f;
	m_fitEnergy=0.0f;
	m_xposLayer.Reset();
	m_yposLayer.Reset();
	m_eneLayer.Reset();
    Create();
}

CalCluster::~CalCluster() {
    Clean();
}

void CalCluster::Clean() {
    m_id = 0;
    m_Esum = 0.0f;
    m_Ecorr = 0.0f;
	m_CsiAlpha=0.0f;
	m_CsiLambda=0.0f;
	m_start=0.0f;
	m_ProfChisq=0.0f;
	m_fitEnergy=0.0f;
	m_xposLayer.Reset();
	m_yposLayer.Reset();
	m_eneLayer.Reset();

    if (m_position) {
        delete m_position;
        m_position = 0;
    }
    
    if (m_direction) {
        delete m_direction;
        m_direction = 0;
    }

    // Don't delete the logs, just
    // clear out the TObjArray,
    // since these logs are probably in 
    // the "master" log list pointed to
    // by CalRecon
    if (m_calLogs) {
        m_calLogs->Clear();
    }
}

void CalCluster::Create() {
    m_calLogs = new TObjArray();
}
