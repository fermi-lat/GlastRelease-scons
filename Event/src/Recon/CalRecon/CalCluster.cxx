
#include "Event/Recon/CalRecon/CalCluster.h"

using namespace Event;

//----------------- CalCluster ------------------

//################################################
CalCluster::CalCluster(double e,Point p)
//################################################
{ 
	int nl = 8;
	m_eneLayer.resize(nl);	
	m_pLayer.resize(nl);	
	ini();
	m_energySum = e;
	m_energyCorrected = m_energySum;
	m_position = p;
}
//################################################
void CalCluster::writeOut() const
//################################################
{
#if 0 // THB: enable this after it is modified to write to the log object, and the values are all defined
	std::cout << "Energy " << m_energySum << " Corrected " << m_energyCorrected;
	std::cout << " " << getPosition().x() << " " << getPosition().y() << " " << getPosition().z();
	std::cout << " " << getDirection().x() << " " << getDirection().y() << " " << getDirection().z();
	std::cout<<"\n";
#endif
}
//------------- private --------------------------
//################################################
void CalCluster::ini()
//################################################
{
	m_energySum       = 0.;
	m_energyCorrected = 0.;

	m_position = Point(0.,0.,0.);
	m_direction = Vector(0.,0.,0);
	int nLayers = m_eneLayer.size();
	for(int i = 0; i<nLayers; i++){
		m_eneLayer[i]=0.;
		m_pLayer[i]=Vector(0.,0.,0.);
	}
}
//----------------- CsIClusterList -----------------

//################################################
void CalClusterCol::clear()
//################################################
{
	int nClusters = num();
	for (int icl = 0; icl < nClusters; icl++) {
		delete m_calClusterCol[icl];
	}
	m_calClusterCol.clear();
}

//------------ private ---------------------------
//################################################
void CalClusterCol::ini()
//################################################
{
	m_calClusterCol.clear();
}
//################################################
void CalClusterCol::writeOut() const
//################################################
{
#if 0 // fix this to write to the log file for debug purposes
    if (m_calClusterCol.size()<=0) return;
	
	std::cout << " --- CalClusterCol  --- " << m_calClusterCol.size() <<"\n";
	for (int i = 0; i < m_calClusterCol.size();i++) {
		m_calClusterCol[i]->writeOut();
	}
#endif
}
