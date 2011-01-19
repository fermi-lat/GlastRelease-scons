
/**
   @file CalXtalsParams.cxx

   @brief Implementation of the CalXtalsParams class.
      
   @author Luca Baldini (luca.baldini@pi.infn.it)
   
   $Revision$
   $Date$
   $Header$
*/


#include "Event/Recon/CalRecon/CalXtalsParams.h"



Event::CalXtalsParams::CalXtalsParams(int numXtals, int numTruncXtals, int numSaturatedXtals,
				      double xtalRawEneSum, double xtalCorrEneSum, double xtalEneMax,
				      double xtalEneRms, double xtalEneSkewness,
				      const Point& centroid) :
  m_numXtals(numXtals),
  m_numTruncXtals(numTruncXtals),
  m_numSaturatedXtals(numSaturatedXtals),
  m_xtalRawEneSum(xtalRawEneSum),
  m_xtalCorrEneSum(xtalCorrEneSum),
  m_xtalEneMax(xtalEneMax),
  m_xtalEneRms(xtalEneRms),
  m_xtalEneSkewness(xtalEneSkewness),
  m_centroid(centroid)
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
  m_xtalEneMax        = 0.;
  m_xtalEneRms        = -1.;
  // Find a sensible value, here.
  m_xtalEneSkewness   = 0.;
  m_centroid          = Point(0.,0.,0.);
}

std::ostream& Event::CalXtalsParams::fillStream( std::ostream& s ) const
{
  s <<
    "Number of xtals = " << m_numXtals << "\n" <<
    "Truncated number of xtals = " << m_numTruncXtals << "\n" <<
    "Number of saturated xtals = " << m_numSaturatedXtals << "\n" <<
    "Raw sum of xtal energies = " << m_xtalRawEneSum << " MeV\n" <<
    "Corrected sum of xtal energies = " << m_xtalCorrEneSum << " MeV\n" <<
    "Maximum xtal energy = " << m_xtalEneMax << " MeV\n" <<
    "Rms of xtal energy distribution = " << m_xtalEneRms << " MeV\n" <<
    "Skewness of xtal energy distribution = " << m_xtalEneSkewness << "\n"
    "Centroid = (" << m_centroid.x() << ", " << m_centroid.y() << ", "
			 << m_centroid.z() << ") mm";

  return s;
}

