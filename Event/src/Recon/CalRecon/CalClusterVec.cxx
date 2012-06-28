

// Implementation file of CalClusterVec class.
//  
// Authors: Luca Baldini.


#include "Event/Recon/CalRecon/CalClusterVec.h"
#include <stdexcept>


Event::CalClusterVec::CalClusterVec()
{
  // Initialize the pointers to the highest-energy and highest gam
  // prob clusters to NULL.
  m_highestEnergyCluster = NULL;
  m_highestGamProbCluster = NULL;
}


void Event::CalClusterVec::push_back(CalCluster* calCluster)
{
  // Keep track of the cluster with the highest energy.
  if (m_highestEnergyCluster == NULL ||
      calCluster->getXtalsParams().getXtalRawEneSum() >
      m_highestEnergyCluster->getXtalsParams().getXtalRawEneSum())
    {
      m_highestEnergyCluster = calCluster;
    }
  // Keep track of the cluster with the highest gamma probability.
  if (m_highestGamProbCluster == NULL ||
      calCluster->getClassParams().getGamProb() >
      m_highestGamProbCluster->getClassParams().getGamProb())
    {
      m_highestGamProbCluster = calCluster;
    }
  // Finally, call the base class push_back() method.
  std::vector<CalCluster*>::push_back(calCluster);
}


