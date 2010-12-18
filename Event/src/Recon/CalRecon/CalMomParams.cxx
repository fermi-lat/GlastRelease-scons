
// Implementation file of CalMomParams class.
//  
// Authors: Luca Baldini, Johan Bregeon.


#include "Event/Recon/CalRecon/CalMomParams.h"


Event::CalMomParams::CalMomParams(double energy, double eneError,
				  const Point& centroid, const CLHEP::HepMatrix& centroidErr,
				  const Vector& axis, const CLHEP::HepMatrix& axisErr,
				  int numIterations, double transRms, double longRms,
				  double longRmsAsym, double longSkewness,
				  double coreEnergyFrac) :
  CalParams(energy, eneError, centroid, centroidErr, axis, axisErr),
  m_numIterations(numIterations),
  m_transRms(transRms),
  m_longRms(longRms),
  m_longRmsAsym(longRmsAsym),
  m_longSkewness(longSkewness),
  m_coreEnergyFrac(coreEnergyFrac)
{
  // Nothing to do, here.
}

Event::CalMomParams::CalMomParams(double energy, double eneError,
				  double xCntrd, double yCntrd, double zCntrd,
				  double cntdxx, double cntdxy, double cntdxz,
				  double cntdyy, double cntdyz, double cntdzz,
				  double xAxis,  double yAxis,  double zAxis,
				  double axsdxx, double axsdxy, double axsdxz,
				  double axsdyy, double axsdyz, double axsdzz,
				  int numIterations, double transRms, double longRms,
				  double longRmsAsym, double longSkewness,
				  double coreEnergyFrac) :
  CalParams(energy, eneError, xCntrd, yCntrd, zCntrd, cntdxx, cntdxy, cntdxz,
	    cntdyy, cntdyz, cntdzz, xAxis, yAxis, zAxis, axsdxx, axsdxy, axsdxz,
	    axsdyy, axsdyz, axsdzz),
  m_numIterations(numIterations),
  m_transRms(transRms),
  m_longRms(longRms),
  m_longRmsAsym(longRmsAsym),
  m_longSkewness(longSkewness),
  m_coreEnergyFrac(coreEnergyFrac)
{
  // Nothing to do, here.
} 

void Event::CalMomParams::clear()
{
  // Call the base class clear method...
  Event::CalParams::clear();

  // ... then reset the additional class members.
  m_numIterations  = 0;
  m_transRms       = -1.;
  m_longRms        = -1.;
  m_longRmsAsym    = -1.;
  m_longSkewness   = -9999.;
  m_coreEnergyFrac = -1.;
}

double Event::CalMomParams::getElongation() const
{
  if (m_transRms > 0.){
    return m_longRms/m_transRms;
  }
  return -1.;
}

std::ostream& Event::CalMomParams::fillStream(std::ostream& s) const
{
  // Call the base class method...
  Event::CalParams::fillStream(s) << "\n";

  // ... then print the additional stuff.
  s <<
    "Number of iterations = " << m_numIterations << "\n" << 
    "Transverse RMS = " << m_transRms << " mm\n" <<
    "Longitudinal RMS = " << m_longRms << " mm\n" <<
    "Longitudinal RMS asymmetry = " << m_longRmsAsym << "\n" <<
    "Longitudinal skewness = " << m_longSkewness << "\n" <<
    "Elongation = " << getElongation() << "\n" <<
    "Core energy fraction = " << m_coreEnergyFrac;  

  return s; 
}
