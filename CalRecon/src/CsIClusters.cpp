
#include "CalRecon/CsIClusters.h"
// #include "Event/messageManager.h"

//----------------- CsICluster ------------------

//################################################
CsICluster::CsICluster(double e,Point p)
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
void CsICluster::writeOut() const
//################################################
{

	std::cout << "Energy " << m_energySum << " Corrected " << m_energyCorrected;
	std::cout << " " << position().x() << " " << position().y() << " " << position().z();
	std::cout << " " << direction().x() << " " << direction().y() << " " << direction().z();
	std::cout<<"\n";
}
//------------- private --------------------------
//################################################
void CsICluster::ini()
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
void CsIClusterList::clear()
//################################################
{
	int nClusters = num();
	for (int icl = 0; icl < nClusters; icl++) {
		delete m_CsIClustersList[icl];
	}
	m_CsIClustersList.clear();
}

//------------ private ---------------------------
//################################################
void CsIClusterList::ini()
//################################################
{
	m_CsIClustersList.clear();
}
//################################################
void CsIClusterList::writeOut() const
//################################################
{
	if (m_CsIClustersList.size()<=0) return;
	
	std::cout << " --- CsIClusterList  --- " << m_CsIClustersList.size() <<"\n";
	for (int i = 0; i < m_CsIClustersList.size();i++) {
		m_CsIClustersList[i]->writeOut();
	}
}
