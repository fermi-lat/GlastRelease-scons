// File and Version information:
// $Header$
//
//  Implementation file of CalCluster and CalClusterCol classes
//  
// Authors:
//
//    Tracy Usher
//
//

#include "Event/Recon/CalRecon/CalParams.h"



Event::CalParams::CalParams(double energy, double eneError,
                            const Point&  centroid, const CLHEP::HepMatrix& centroidErr,
                            const Vector& axis,     const CLHEP::HepMatrix& axisErr) :
  m_energy(energy),
  m_eneError(eneError),
  m_clusterCentroid(centroid),
  m_clusterAxis(axis)
{
  setCentroidErrs(centroidErr);
  setAxisErrs(axisErr);
}

Event::CalParams::CalParams(double energy, double eneError,
                            double xCntrd, double yCntrd, double zCntrd,
                            double cntdxx, double cntdxy, double cntdxz,
			    double cntdyy, double cntdyz, double cntdzz,
                            double xAxis,  double yAxis,  double zAxis,
                            double axsdxx, double axsdxy, double axsdxz,
			    double axsdyy, double axsdyz, double axsdzz) :
  m_energy(energy),
  m_eneError(eneError),
  m_cenxx(cntdxx),
  m_cenxy(cntdxy),
  m_cenxz(cntdxz),
  m_cenyy(cntdyy),
  m_cenyz(cntdyz),
  m_cenzz(cntdzz),
  m_axisxx(axsdxx),
  m_axisxy(axsdxy),
  m_axisxz(axsdxz),
  m_axisyy(axsdyy),
  m_axisyz(axsdyz),
  m_axiszz(axsdzz)
{
  m_clusterCentroid = Point(xCntrd,yCntrd,zCntrd);
  m_clusterAxis     = Vector(xAxis,yAxis,zAxis);
}

Event::CalParams::CalParams(double energy, double eneError,
                            const Point&  centroid, const CLHEP::HepMatrix& centroidErr) :
  m_energy(energy),
  m_eneError(eneError),
  m_clusterCentroid(centroid)
{
  setCentroidErrs(centroidErr);
  clearAxis();
}

Event::CalParams::CalParams(double energy, double eneError,
                            double xCntrd, double yCntrd, double zCntrd,
                            double cntdxx, double cntdxy, double cntdxz,
			    double cntdyy, double cntdyz, double cntdzz) :
  m_energy(energy),
  m_eneError(eneError),
  m_cenxx(cntdxx),
  m_cenxy(cntdxy),
  m_cenxz(cntdxz),
  m_cenyy(cntdyy),
  m_cenyz(cntdyz),
  m_cenzz(cntdzz)
{
  m_clusterCentroid = Point(xCntrd,yCntrd,zCntrd);
  clearAxis();
}

void Event::CalParams::clearEnergy()
{
  m_energy   = 0.;
  m_eneError = 0.;
}

void Event::CalParams::clearCentroid()
{
  m_clusterCentroid = Point(0.,0.,0.);

  m_cenxx = 0.;
  m_cenxy = 0.;
  m_cenxz = 0.;
  m_cenyy = 0.;
  m_cenyz = 0.;
  m_cenzz = 0.;
}

void Event::CalParams::clearAxis()
{
  m_clusterAxis = Vector(0.,0.,0.);

  m_axisxx = 0.;
  m_axisxy = 0.;
  m_axisxz = 0.;
  m_axisyy = 0.;
  m_axisyz = 0.;
  m_axiszz = 0.;
}

void Event::CalParams::clear()
{
  clearEnergy();
  clearCentroid();
  clearAxis();
}

void Event::CalParams::setCentroidErrs(const CLHEP::HepMatrix& errs)
{
  m_cenxx = errs(1,1);
  m_cenxy = errs(1,2);
  m_cenxz = errs(1,3);
  m_cenyy = errs(2,2);
  m_cenyz = errs(2,3);
  m_cenzz = errs(3,3);
}

void Event::CalParams::setAxisErrs(const CLHEP::HepMatrix& errs)
{
  m_axisxx = errs(1,1);
  m_axisxy = errs(1,2);
  m_axisxz = errs(1,3);
  m_axisyy = errs(2,2);
  m_axisyz = errs(2,3);
  m_axiszz = errs(3,3);
}

CLHEP::HepMatrix Event::CalParams::getCentroidErrs() const
{
  CLHEP::HepMatrix errs(3,3,0);
  
  errs(1,1) = m_cenxx;
  errs(1,2) = errs(2,1) = m_cenxy;
  errs(1,3) = errs(3,1) = m_cenxz;
  errs(2,2) = m_cenyy;
  errs(2,3) = errs(3,2) = m_cenyz;
  errs(3,3) = m_cenzz;
  
  return errs;
}

CLHEP::HepMatrix Event::CalParams::getAxisErrs() const
{
  CLHEP::HepMatrix errs(3,3,0);
  
  errs(1,1) = m_axisxx;
  errs(1,2) = errs(2,1) = m_axisxy;
  errs(1,3) = errs(3,1) = m_axisxz;
  errs(2,2) = m_axisyy;
  errs(2,3) = errs(3,2) = m_axisyz;
  errs(3,3) = m_axiszz;
  
  return errs;
}

std::ostream& Event::CalParams::fillStream( std::ostream& s ) const
{
  s <<
    "Energy = " << m_energy << " +- " << m_eneError << " MeV\n" <<
    "Centroid = (" << m_clusterCentroid.x() << ", " << m_clusterCentroid.y() << ", "
		<< m_clusterCentroid.z() << ") mm\n" <<
    "Centroid covariance matrix:\n" <<
    "| " << m_cenxx  << "  " << m_cenxy << "  " << m_cenxz << " |\n" <<
    "| " << m_cenxy  << "  " << m_cenyy << "  " << m_cenyz << " |\n" <<
    "| " << m_cenxz  << "  " << m_cenyz << "  " << m_cenzz << " |\n" <<
    "Axis = (" << m_clusterAxis.x() << ", " << m_clusterAxis.y() << ", "
		  << m_clusterAxis.z() << ")\n" <<
    "Axis covariance matrix:\n" <<
    "| " << m_axisxx  << "  " << m_axisxy << "  " << m_axisxz << " |\n" <<
    "| " << m_axisxy  << "  " << m_axisyy << "  " << m_axisyz << " |\n" <<
    "| " << m_axisxz  << "  " << m_axisyz << "  " << m_axiszz << " |";

    return s; 
}

