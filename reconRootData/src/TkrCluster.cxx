#include "reconRootData/TkrCluster.h"
#include <iostream>


ClassImp(TkrCluster)

TkrCluster::TkrCluster() {
    Clear();
}

TkrCluster::TkrCluster(UInt_t id, UInt_t iplane, TkrCluster::view v,
                       UInt_t istrip0, UInt_t istripf, const TVector3& position, 
					   Double_t ToT, UInt_t flag, UInt_t tower)
{
	// Purpose and method: makes a cluster with attributes
	// Input:  id is the sequential cluster number
	//         v  is the measured view (0->X, 1->Y)
	//         iplane is the bilayer number (0 at the front)
	//         istrip0 is the first strip in the cluster
	//         istripf is the last strip in the cluster
	//         ToT is the time-over-threshold for this cluster
	//         flag is
	//         tower is the tower number

	Clear();
	
    m_id     = id;
	m_view   = v;
	
	m_plane  = iplane;
	
	m_strip0 = istrip0;
	m_stripf = istripf;
	
	m_position = position;
	
	m_ToT    = ToT;
    m_tower  = tower;
	
    m_flag   = flag;
}


TkrCluster::~TkrCluster() {
    Clear();
}



void TkrCluster::Clear(Option_t *option) {
	m_id = 0;
	m_tower = 0;
	m_plane=0;
    m_flag = 0;
	m_strip0 = 0;
    m_stripf = 0;
    m_view = TkrCluster::X;
    m_ToT = 0.0;
    m_position = TVector3(0.0, 0.0, 0.0);
}

void TkrCluster::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;

	cout << "(Id, Tower, Plane): (" << m_id << ", " << m_tower 
		<< ", " << m_plane << ")" << endl;
	cout << "First Strip, Last Strip " << m_strip0 << ", " << m_stripf << endl;
	cout << "ToT, view, position(x, y, z) " << m_ToT << " " << m_view << " ("
		<< m_position.X() << ", " << m_position.Y() << ", " << m_position.Z() 
		<< ")" << endl;
}
