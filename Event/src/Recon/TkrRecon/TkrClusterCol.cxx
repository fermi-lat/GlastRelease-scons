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

using namespace Event;

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

int TkrClusterCol::nHits(TkrCluster::view v, int iplane) const
{
    if ((v == TkrCluster::X || v == TkrCluster::Y) && iplane >= 0 && iplane < NPLANES)
    {
        return (int) getHits(v,iplane).size();
    }
    else 
    {
        return 0;
    }
}

TkrCluster* TkrClusterCol::getHit(int i) const 
{
    // this used to just return the ith cluster...
    // now it really checks the id's, but tries the ith cluster first 
    //   (which will be the one, unless someone sorts the list of pointers)
    TkrCluster* ptr = 0;
    if (i==m_clustersList[i]->id()) {
        ptr = m_clustersList[i];
    } else {
        for (int j = 0; j< m_clustersList.size(); j++) {
            if (m_clustersList[j]->id()==i) ptr = m_clustersList[j];
        }
    }
    return ptr;
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


