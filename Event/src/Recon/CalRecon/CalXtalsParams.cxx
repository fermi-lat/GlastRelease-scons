
//
//  Implementation file of the CalXtalsParams container.
//  
// Authors:
//
//    Tracy Usher
//
//

#include "Event/Recon/CalRecon/CalXtalsParams.h"



Event::CalXtalsParams::CalXtalsParams(int numXtals, int numTruncXtals, int numSaturatedXtals,
				      double xtalRawEneSum, double xtalCorrEneSum,
				      double xtalEneRms, double xtalEneSkewness) :
  m_numXtals(numXtals),
  m_numTruncXtals(numTruncXtals),
  m_numSaturatedXtals(numSaturatedXtals),
  m_xtalRawEneSum(xtalRawEneSum),
  m_xtalCorrEneSum(xtalCorrEneSum),
  m_xtalEneRms(xtalEneRms),
  m_xtalEneSkewness(xtalEneSkewness)
{
  // Nothing to do, here.
}

Event::CalXtalsParams::CalXtalsParams(int numTruncXtals, int numSaturatedXtals)
{
  clear();
  m_numTruncXtals     = numTruncXtals;
  m_numSaturatedXtals = numSaturatedXtals;
}

void Event::CalXtalsParams::clear()
{
  m_numXtals          = 0;
  m_numTruncXtals     = 0;
  m_numSaturatedXtals = 0;
  m_xtalRawEneSum     = 0.;
  m_xtalCorrEneSum    = 0.;
  m_xtalEneRms        = -1.;
  // Find a sensible value, here.
  m_xtalEneSkewness   = 0.;
}

std::ostream& Event::CalXtalsParams::fillStream( std::ostream& s ) const
{
  s <<
    "Number of xtals = " << m_numXtals << "\n" <<
    "Truncated number of xtals = " << m_numTruncXtals << "\n" <<
    "Number of saturated xtals = " << m_numSaturatedXtals << "\n" <<
    "Raw sum of xtal energies = " << m_xtalRawEneSum << "\n" <<
    "Corrected sum of xtal energies = " << m_xtalCorrEneSum << "\n" <<
    "Rms of xtal energy distribution = " << m_xtalEneRms << "\n" <<
    "Skewness of xtal energy distribution = " << m_xtalEneSkewness;
  
  return s; 
}

