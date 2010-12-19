// File and Version information:
// $Header$
//
//  Implementation file of CalCluster and CalClusterCol classes
//  
// Authors:
//
//    Luca Baldini
//    Johan Bregeon
//

#include "Event/Recon/CalRecon/CalCluster.h"

void Event::CalCluster::initialize(const CalMSTreeParams& mstParams,
				   const CalFitParams& fitParams,
				   const CalMomParams& momParams,
				   const std::map <std::string, double>& classProb,
				   int numSaturatedXtals, int numTruncXtals)
{
  m_mstParams         = mstParams;
  m_fitParams         = fitParams;
  m_momParams         = momParams;
  m_classesProb       = classProb;
  m_numSaturatedXtals = numSaturatedXtals;
  m_numTruncXtals     = numTruncXtals;
}

/// TBD: rename to clear()
void Event::CalCluster::iniCluster()
{
  m_mstParams          = CalMSTreeParams();
  m_fitParams          = CalFitParams();
  m_momParams          = CalMomParams();
  // Create the map and initialize gam prob to -1 -- needed ? TBD
  m_classesProb        = std::map <std::string, double>();
  m_classesProb["gam"] = -1.;
  m_statusBits         = 0;
  m_numSaturatedXtals  = 0;
  m_numTruncXtals      = 0; 
}

double Event::CalCluster::getTopologyProb(std::string top) const
{
  if(m_classesProb.count(top))
    return m_classesProb.find(top)->second;
  else
    return -1;
}

void Event::CalCluster::writeOut(MsgStream& stream) const
{
  stream << "Energy " << m_momParams.getEnergy();
  stream << " No.Trunc Xtals " << m_numTruncXtals;
  stream << " Gam prob " << getGamProb();
  stream << " " << getPosition().x() 
	 << " " << getPosition().y() 
	 << " " << getPosition().z();
  stream << " " << getDirection().x() 
	 << " " << getDirection().y() 
	 << " " << getDirection().z();
  stream << endreq;
}

std::ostream& Event::CalCluster::fillStream(std::ostream& s) const
{
  s <<
    "--- Moments analysis output ---\n" << m_momParams << "\n";
  return s; 
}
