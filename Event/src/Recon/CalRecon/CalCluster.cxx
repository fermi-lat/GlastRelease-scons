
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
void CalCluster::writeOut(MsgStream& stream) const
//################################################
{
#if 1 // THB: enable this after it is modified to write to the log object, and the values are all defined
	stream << "Energy " << m_energySum << " Corrected " << m_energyCorrected;
	stream << " " << getPosition().x() << " " << getPosition().y() << " " << getPosition().z();
	stream << " " << getDirection().x() << " " << getDirection().y() << " " << getDirection().z();
	stream << endreq;
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
void CalClusterCol::writeOut(MsgStream& stream) const
//################################################
{
#if 1 // fix this to write to the log file for debug purposes
    if (m_calClusterCol.size()<=0) return;
	
	stream << " --- CalClusterCol  --- " << m_calClusterCol.size() << endreq;
	for (int i = 0; i < m_calClusterCol.size();i++) {
		m_calClusterCol[i]->writeOut(stream);
	}
#endif
}
