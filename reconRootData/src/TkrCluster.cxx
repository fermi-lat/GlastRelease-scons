#include "reconRootData/TkrCluster.h"
#include <iostream>


ClassImp(TkrCluster)

TkrCluster::TkrCluster() {
    Clear();
}

TkrCluster::TkrCluster(commonRootData::TkrId tkrId, UInt_t istrip0, UInt_t istripf, 
                       const TVector3& position, Float_t ToT, UInt_t status, UInt_t nBad)
{
	// Purpose and method: makes a cluster with attributes
	// Input:  tkrId is a valid TkrId which identifies where the cluster is
    //         strip0 is the first strip in the cluster
    //         stripf is the last strip in the cluster
    //         position is the position of the cluster
    //         ToT is the time-over-threshold for this cluster
	//         flag tells if the cluster is "in use" or not
    //         id is the sequential cluster number

	Clear();
	
    m_tkrId    = tkrId;
    m_strip0   = istrip0;
    m_stripf   = istripf;	
    m_position = position;	
    m_ToT      = ToT;	
    m_status   = status;
    m_nBad     = nBad;
}


TkrCluster::~TkrCluster() 
{
    Clear();
}



void TkrCluster::Clear(Option_t* /* option */) 
{
    m_tkrId    = commonRootData::TkrId();
    m_status   = 0;
    m_strip0   = 0;
    m_stripf   = 0;
    m_ToT      = 0.0;
    m_position = TVector3(0.0, 0.0, 0.0);
    m_nBad     = 0;
}

void TkrCluster::Print(Option_t *option) const 
{
    TObject::Print(option);
    using namespace std;

	cout << "(TkrId): (" << m_tkrId << endl;
	cout << "First Strip, Last Strip " << m_strip0 << ", " << m_stripf << endl;
	cout << "ToT, position(x, y, z) " << m_ToT  << " ("
		<< m_position.X() << ", " << m_position.Y() << ", " << m_position.Z() 
		<< ")" << endl;
}
