
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
                                  double fullLength, double dEdxSpread, double minGhostDoca) :
  CalParams(energy, eneError, centroid, centroidErr, axis, axisErr),
  m_numIterations(numIterations),
  m_numCoreXtals(numCoreXtals),
  m_numXtals(numXtals),
  m_transRms(transRms),
  m_longRms(longRms),
  m_longRmsAsym(longRmsAsym),
  m_longSkewness(longSkewness),
  m_coreEnergyFrac(coreEnergyFrac),
  m_fullLength(fullLength),
  m_dEdxSpread(dEdxSpread),
  m_minGhostDoca(minGhostDoca)
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
                                  double fullLength, double dEdxSpread, double minGhostDoca) :
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
  m_fullLength(fullLength),
  m_dEdxSpread(dEdxSpread),
  m_minGhostDoca(minGhostDoca)
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

Event::CalMomParams::CalMomParams(double energy, double eneError,
                                  double xCntrd, double yCntrd, double zCntrd,
                                  double xAxis,  double yAxis,  double zAxis)
{
  clear();
  setEnergy(energy);
  setEnergyErr(eneError);
  setCentroid( Point(xCntrd,yCntrd,zCntrd) );
  setAxis( Vector(xAxis,yAxis,zAxis) );
}

Event::CalMomParams::CalMomParams(double energy, double eneError,
                                  const Point& centroid, const Vector& axis)
{
  clear();
  setEnergy(energy);
  setEnergyErr(eneError);
  setCentroid(centroid);
  setAxis(axis);
}

Event::CalMomParams::CalMomParams(double energy, double eneError,
                                  double xCntrd, double yCntrd, double zCntrd)
{
  clear();
  setEnergy(energy);
  setEnergyErr(eneError);
  setCentroid( Point(xCntrd,yCntrd,zCntrd) );
}

Event::CalMomParams::CalMomParams(double energy, double eneError,
                                  const Point& centroid)
{
  clear();
  setEnergy(energy);
  setEnergyErr(eneError);
  setCentroid(centroid);
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
  m_fullLength     = 0.;
  m_dEdxSpread     = 0.;
  m_minGhostDoca   = -1.;
}

double Event::CalMomParams::getElongation() const
{
  if ( getTransRms() > 0. ) {
    return getLongRms()/getTransRms();
  }
  return -1.;
}

double Event::CalMomParams::getdEdxAverage() const
{
  if ( getFullLength() > 0. ) {
    return getEnergy()/getFullLength();
  }
  return -1.;
}

std::ostream& Event::CalMomParams::fillStream(std::ostream& s) const
{
  // Call the base class method...
  Event::CalParams::fillStream(s) << "\n";

  // ... then print the additional stuff.
  s <<
    "Number of iterations = " << getNumIterations() << "\n" << 
    getNumCoreXtals() << "/" << getNumXtals() <<
    " xtal(s) used for the final centroid/axis\n" << 
    "Transverse RMS = " << getTransRms() << " mm\n" <<
    "Longitudinal RMS = " << getLongRms() << " mm\n" <<
    "Longitudinal RMS asymmetry = " << getLongRmsAsym() << "\n" <<
    "Longitudinal skewness = " << getLongSkewness() << "\n" <<
    "Elongation = " << getElongation() << "\n" <<
    "Core energy fraction = " << getCoreEnergyFrac() << "\n" <<
    "Cluster full length = " << getFullLength() << " X0\n" <<
    "Average dE/dx = " << getdEdxAverage() << " MeV/X0 (+- " << getdEdxSpread() << ")\n" <<
    "Minimum DOCA to Ghost Track = " << getMinGhostDoca() << " mm";

  return s; 
}
