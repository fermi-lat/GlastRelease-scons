#include "reconRootData/CalCluster.h"


ClassImp(CalCluster)


//################################################
CalCluster::CalCluster(Double_t e,TVector3 p)
//################################################
{ 
	int nl = 8;
	m_eneLayer.resize(nl);	
	m_pLayer.resize(nl);	
	ini();
	m_energySum = e;
	m_energyCorrected = m_energySum;
	m_position = p;
};
//################################################
void CalCluster::writeOut() const
//################################################
{
#if 0 // THB: enable this after it is modified to write to the log object, and the values are all defined
	std::cout << "Energy " << m_energySum << " Corrected " << m_energyCorrected;
	std::cout << " " << position().x() << " " << position().y() << " " << position().z();
	std::cout << " " << direction().x() << " " << direction().y() << " " << direction().z();
	std::cout<<"\n";
#endif
};
//------------- private --------------------------
//################################################
void CalCluster::ini()
//################################################
{
	m_energySum       = 0.;
	m_energyCorrected = 0.;

	m_position = TVector3(0.,0.,0.);
	m_direction = TVector3(0.,0.,0);
	int nLayers = m_eneLayer.size();
	for(int i = 0; i<nLayers; i++){
		m_eneLayer[i]=0.;
		m_pLayer[i]=TVector3(0.,0.,0.);
	}
};
