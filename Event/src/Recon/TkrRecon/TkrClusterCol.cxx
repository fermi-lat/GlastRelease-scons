//      $Header$
//
// Description:
//      TkrClusterCol is a container for Tkr clusters, and has the methods
//      for accessing the cluster information.
//
// Author(s):
//      Tracy Usher     
//      Leon Rochester     



#include "Event/Recon/TkrRecon/TkrClusterCol.h"

using namespace TkrRecon;

TkrClusterCol::TkrClusterCol()
{
    //Initialize the cluster lists...
    ini();
}


TkrClusterCol::~TkrClusterCol()
{
    // This deletes all the clusters
	clear();
	
    return;
}

void TkrClusterCol::addCluster(TkrCluster* cl)
{
    // Purpose and Method: Adds a cluster to the cluster lists
    // Inputs:  cl is the cluster to be added
	m_clustersList.push_back(cl);
	int iview = TkrCluster::viewToInt(cl->v());
	m_clustersByPlaneList[iview][cl->plane()].push_back(cl);
}

void TkrClusterCol::clear()
{
    // Purpose and Method: deletes the clusters
	int nhits = m_clustersList.size();
	for (int ihit = 0; ihit < nhits; ihit++) {
		delete m_clustersList[ihit];
	}
	ini();
}
void TkrClusterCol::ini()
{
    // Purpose and Method: clears all the cluster lists
    // Inputs:  None
	
    // this "clear" is the clear method of std::vector
    //   not TkrClusterCol::clear!
    m_clustersList.clear();
	for (int iview = 0; iview < NVIEWS; iview++) {
		for (int iplane = 0; iplane < NPLANES; iplane++) {
			m_clustersByPlaneList[iview][iplane].clear();
		}
	}
}


void TkrClusterCol::writeOut(MsgStream& log) const
{
	// Purpose: writes out the information about the clusters
	// Method: calls writeOut() for each cluster.
	if (nHits()<=0) return;
	
	for (int ihit = 0; ihit < nHits(); ihit++) {
		m_clustersList[ihit]->writeOut(log);
	}
}


