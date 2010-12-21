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
  std::stringstream s;
  fillStream(s);
  stream << s;
}

std::ostream& Event::CalCluster::fillStream(std::ostream& s) const
{
  s <<
    "Cal cluster status bits: " << m_statusBits << "\n" <<
    "Number of saturated xtals: " << m_numSaturatedXtals << "\n" <<
    "Truncated number of xtals: " << m_numTruncXtals << "\n" <<
    "----------------------------------------------------\n" <<
    "---------- Output from the MST clustering ----------\n" << 
    "----------------------------------------------------\n" << m_mstParams << "\n" <<
    "----------------------------------------------------\n" <<
    "----------- Output from the fitting tool -----------\n" << 
    "----------------------------------------------------\n" << m_fitParams << "\n" <<
    "----------------------------------------------------\n" <<
    "--------- Output from the moments analysis ---------\n" << 
    "----------------------------------------------------\n" << m_momParams << "\n" <<
    "------------------------------------------------------";
  return s; 
}
