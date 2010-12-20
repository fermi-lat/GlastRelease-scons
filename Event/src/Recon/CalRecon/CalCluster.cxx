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
				   const CalClassParams& classParams,
				   int numSaturatedXtals, int numTruncXtals)
{
  m_mstParams         = mstParams;
  m_fitParams         = fitParams;
  m_momParams         = momParams;
  m_classParams       = classParams;
  m_numSaturatedXtals = numSaturatedXtals;
  m_numTruncXtals     = numTruncXtals;
}

/// TBD: rename to clear()
void Event::CalCluster::iniCluster()
{
  m_mstParams         = CalMSTreeParams();
  m_fitParams         = CalFitParams();
  m_momParams         = CalMomParams();
  m_classParams       = CalClassParams();
  m_statusBits        = 0;
  m_numSaturatedXtals = 0;
  m_numTruncXtals     = 0; 
}

double Event::CalCluster::getClassProb(const std::string& className) const
{
  return m_classParams.getProb(className);
}

double Event::CalCluster::getGamProb() const
{
  return getClassProb("gam");
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
