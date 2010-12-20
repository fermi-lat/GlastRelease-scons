
// Implementation file of CalMomParams class.
//  
// Authors: Luca Baldini, Johan Bregeon.


#include "Event/Recon/CalRecon/CalMomParams.h"


Event::CalMomParams::CalMomParams(double energy, double eneError,
				  const Point& centroid, const CLHEP::HepMatrix& centroidErr,
				  const Vector& axis, const CLHEP::HepMatrix& axisErr,
				  int numIterations, int numCoreXtals, int numXtals,
				  double transRms, double longRms, double longRmsAsym,
				  double longSkewness, double coreEnergyFrac,
				  double dEdxAverage, double dEdxSpread) :
  CalParams(energy, eneError, centroid, centroidErr, axis, axisErr),
  m_numIterations(numIterations),
  m_numCoreXtals(numCoreXtals),
  m_numXtals(numXtals),
  m_transRms(transRms),
  m_longRms(longRms),
  m_longRmsAsym(longRmsAsym),
  m_longSkewness(longSkewness),
  m_coreEnergyFrac(coreEnergyFrac),
  m_dEdxAverage(dEdxAverage),
  m_dEdxSpread(dEdxSpread)
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
				  int numIterations, int numCoreXtals, int numXtals,
				  double transRms, double longRms, double longRmsAsym,
				  double longSkewness, double coreEnergyFrac,
				  double dEdxAverage, double dEdxSpread) :
  CalParams(energy, eneError, xCntrd, yCntrd, zCntrd, cntdxx, cntdxy, cntdxz,
	    cntdyy, cntdyz, cntdzz, xAxis, yAxis, zAxis, axsdxx, axsdxy, axsdxz,
	    axsdyy, axsdyz, axsdzz),
  m_numIterations(numIterations),
  m_numCoreXtals(numCoreXtals),
  m_numXtals(numXtals),
  m_transRms(transRms),
  m_longRms(longRms),
  m_longRmsAsym(longRmsAsym),
  m_longSkewness(longSkewness),
  m_coreEnergyFrac(coreEnergyFrac),
  m_dEdxAverage(dEdxAverage),
  m_dEdxSpread(dEdxSpread)
{
  // Nothing to do, here.
}

Event::CalMomParams::CalMomParams(double energy, double eneError,
				  double xCntrd, double yCntrd, double zCntrd,
				  double cntdxx, double cntdxy, double cntdxz,
				  double cntdyy, double cntdyz, double cntdzz,
				  double xAxis,  double yAxis,  double zAxis,
				  double axsdxx, double axsdxy, double axsdxz,
				  double axsdyy, double axsdyz, double axsdzz) :
  CalParams(energy, eneError, xCntrd, yCntrd, zCntrd, cntdxx, cntdxy, cntdxz,
	    cntdyy, cntdyz, cntdzz, xAxis, yAxis, zAxis, axsdxx, axsdxy, axsdxz,
	    axsdyy, axsdyz, axsdzz)
{
  // Initialize the members which are still undefined.
  clearMomParams();
}

void Event::CalMomParams::clear()
{
  // Call the base class clear method...
  Event::CalParams::clear();

  // ... then reset the additional class members.
  clearMomParams();
}

// TBD: define sensible initialization values, here.
// Need to move the moment normalization into the moments analysis before you can put
// negative numbers, here---in order not to get in troubles with sqrt() in AnalysisNtuple
// later on.
// Make sure the values are reflected on the ROOT side.
void Event::CalMomParams::clearMomParams()
{
  m_numIterations  = 0;
  m_numCoreXtals   = 0;
  m_numXtals       = 0;
  m_transRms       = 0.;
  m_longRms        = 0.;
  m_longRmsAsym    = 0.;
  m_longSkewness   = 0.;
  m_coreEnergyFrac = 0.;
  m_dEdxAverage    = 0.;
  m_dEdxSpread     = 0.;
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
    m_numCoreXtals << "/" << m_numXtals << " used for the centroid/axis\n" << 
    "Transverse RMS = " << m_transRms << " mm\n" <<
    "Longitudinal RMS = " << m_longRms << " mm\n" <<
    "Longitudinal RMS asymmetry = " << m_longRmsAsym << "\n" <<
    "Longitudinal skewness = " << m_longSkewness << "\n" <<
    "Elongation = " << getElongation() << "\n" <<
    "Core energy fraction = " << m_coreEnergyFrac << "\n" <<
    "Average dE/dx = " << m_dEdxAverage << " MeV/X0 (+- " << m_dEdxSpread << ")";  

  return s; 
}
