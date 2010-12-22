
// Implementation file of CalFitParams class.
//  
// Authors: Luca Baldini, Johan Bregeon.


#include "Event/Recon/CalRecon/CalFitParams.h"


Event::CalFitParams::CalFitParams(double energy, double eneError,
				  const Point& centroid, const CLHEP::HepMatrix& centroidErr,
				  const Vector& axis, const CLHEP::HepMatrix& axisErr,
				  int numFitLayers, double chiSquare) :
  CalParams(energy, eneError, centroid, centroidErr, axis, axisErr),
  m_nFitLayers(numFitLayers),
  m_chiSquare(chiSquare)
{
  // Nothing to do, here.
}

Event::CalFitParams::CalFitParams(double energy, double eneError,
				  double xCntrd, double yCntrd, double zCntrd,
				  double cntdxx, double cntdxy, double cntdxz,
				  double cntdyy, double cntdyz, double cntdzz,
				  double xAxis,  double yAxis,  double zAxis,
				  double axsdxx, double axsdxy, double axsdxz,
				  double axsdyy, double axsdyz, double axsdzz,
				  int numFitLayers, double chiSquare) :
  CalParams(energy, eneError, xCntrd, yCntrd, zCntrd, cntdxx, cntdxy, cntdxz,
	    cntdyy, cntdyz, cntdzz, xAxis, yAxis, zAxis, axsdxx, axsdxy, axsdxz,
	    axsdyy, axsdyz, axsdzz),
  m_nFitLayers(numFitLayers),
  m_chiSquare(chiSquare)
{
  // Nothing to do, here.
}

Event::CalFitParams::CalFitParams(int numFitLayers, double chiSquare,
				  const Point& centroid, const CLHEP::HepMatrix& centroidErr,
				  const Vector& axis, const CLHEP::HepMatrix& axisErr) :
  CalParams(-1., -1., centroid, centroidErr, axis, axisErr),
  m_nFitLayers(numFitLayers),
  m_chiSquare(chiSquare)
{
  // Nothing to do, here.
}

Event::CalFitParams::CalFitParams(int numFitLayers, double chiSquare,
				  double xCntrd, double yCntrd, double zCntrd,
				  double cntdxx, double cntdxy, double cntdxz,
				  double cntdyy, double cntdyz, double cntdzz,
				  double xAxis,  double yAxis,  double zAxis,
				  double axsdxx, double axsdxy, double axsdxz,
				  double axsdyy, double axsdyz, double axsdzz) :
  CalParams(-1., -1., xCntrd, yCntrd, zCntrd, cntdxx, cntdxy, cntdxz,
	    cntdyy, cntdyz, cntdzz, xAxis, yAxis, zAxis, axsdxx, axsdxy, axsdxz,
	    axsdyy, axsdyz, axsdzz),
  m_nFitLayers(numFitLayers),
  m_chiSquare(chiSquare)
{
  // Nothing to do, here.
}

Event::CalFitParams::CalFitParams(double energy, double eneError,
				  double xCntrd, double yCntrd, double zCntrd,
				  double xAxis,  double yAxis,  double zAxis,
				  int numFitLayers, double chiSquare)
{
  clear();
  setEnergy(energy);
  setEnergyErr(eneError);
  setCentroid( Point(xCntrd,yCntrd,zCntrd) );
  setAxis( Vector(xAxis,yAxis,zAxis) );
  setFitLayers(numFitLayers);
  setChiSquare(chiSquare);
}

Event::CalFitParams::CalFitParams(double energy, double eneError,
				  double xCntrd, double yCntrd, double zCntrd,
				  int numFitLayers, double chiSquare)
{
  clear();
  setEnergy(energy);
  setEnergyErr(eneError);
  setCentroid( Point(xCntrd,yCntrd,zCntrd) );
  setFitLayers(numFitLayers);
  setChiSquare(chiSquare);
}

Event::CalFitParams::CalFitParams(double energy, double eneError,
				  const Point& centroid, int numFitLayers, double chiSquare)
{
  clear();
  setEnergy(energy);
  setEnergyErr(eneError);
  setCentroid(centroid);
  setFitLayers(numFitLayers);
  setChiSquare(chiSquare);
}

void Event::CalFitParams::clear()
{
  // Call the base class clear method...
  Event::CalParams::clear();

  // ... then reset the additional class members.
  clearFitParams();
}

void Event::CalFitParams::clearFitParams()
{
  m_nFitLayers = 0;
  m_chiSquare  = 0.;
}

std::ostream& Event::CalFitParams::fillStream(std::ostream& s) const
{
  // Call the base class method...
  Event::CalParams::fillStream(s) << "\n";

  // ... then print the additional stuff.
  s <<
    "Number of layers for the fit: " << m_nFitLayers << "\n" <<
    "Fit chisquare: " << m_chiSquare;
  
  return s; 
}
