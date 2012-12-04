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

CLHEP::HepMatrix Event::CalParams::getMomErrsAcdRep() const
{ 
  // @brief Return Cal covariance matrix in 5x5 Acd representation
  // Acd rep. is (Vx, Vy, Vz, X0, Y0)
  // this means simply code a matrix with m_axisIJ and m_cenIJ
  // see AcdRecon/doc/formulae.tex for definitions
  
  CLHEP::HepMatrix acdErrs(5,5,0);

  // axis part
  acdErrs(1,1) = m_axisxx;
  acdErrs(1,2) = acdErrs(2,1) = m_axisxy;
  acdErrs(1,3) = acdErrs(3,1) = m_axisxz;
  acdErrs(2,2) = m_axisyy;
  acdErrs(2,3) = acdErrs(3,2) = m_axisyz;
  acdErrs(3,3) = m_axiszz;
  // centroid part
  acdErrs(4,4) = m_cenxx;
  acdErrs(5,5) = m_cenyy;
  acdErrs(4,5) = acdErrs(5,4) = m_cenxy;

  return acdErrs;
}

CLHEP::HepMatrix Event::CalParams::getMomErrsTkrRep() const
{
  // @brief Return Cal covariance matrix in 4x4 Tkr representation
  // Tkr rep. is (X0, Sx, Y0, Sy) with Si = Vi/Vz 
  // The conversion Cov_Tkr = B*Cov_Acd*B.T() with B 4x5
  // see AcdRecon/doc/formulae.tex for definitions

  CLHEP::HepMatrix acdErrs(5,5,0);
  acdErrs = getMomErrsAcdRep();
  CLHEP::HepMatrix tkrErrs(4,4,0);

  CLHEP::HepMatrix B(4,5,0);
  if (m_clusterAxis.z() !=0){
  
    B(1,4) = 1.;
    B(2,1) = 1./m_clusterAxis.z();
    B(2,3) = -1.*m_clusterAxis.x()/(m_clusterAxis.z() * m_clusterAxis.z());
    B(3,5) = 1.;
    B(4,2) = 1./m_clusterAxis.z();
    B(4,3) = -1.*m_clusterAxis.y()/(m_clusterAxis.z() * m_clusterAxis.z());
    tkrErrs = B * acdErrs * B.T();
  }
  return tkrErrs;
}

Point Event::CalParams::getCorCentroid(Vector vaxis)
{
  double rawenergy = m_energy;
  double p[3];
  p[0] = m_clusterCentroid.x();
  p[1] = m_clusterCentroid.y();
  p[2] = m_clusterCentroid.z();
  double v[3];
  v[0] = vaxis.x();
  v[1] = vaxis.y();
  v[2] = vaxis.z();
  if(v[2]>0)
    {
      v[0] = -v[0];
      v[1] = -v[1];
      v[2] = -v[2];
    }
  double mynorm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  double pc[3];
  int i;
  for(i=0;i<3;++i) pc[i] = p[i];
  //
  if(rawenergy<=0 || mynorm==0) return m_clusterCentroid;

  v[0] /= mynorm;
  v[1] /= mynorm;
  v[2] /= mynorm;

  double loge = log10(rawenergy);
  //
  double x,xtow,xtowc,xlog,xrec;
  int itow,ilog;
  //
  double myamp0,myamp,mysig;
  double par0,par1;
  double myxmax,mymax,delta;
  for(i=0;i<2;++i)
    {
      x = p[i];
      if(fabs(x)>1.5*374.5+5*27.84) continue;
      //
      itow = (int)floor((x+749)/374.5);
      xtow = x+749-374.5*(double)itow;
      xtowc = xtow-374.5/2;
      if(fabs(xtowc)>5*27.84) continue;
      //
      ilog = (int)floor((xtow+7.63)/27.84);
      xlog = (xtow+7.63)/27.84-(double)ilog;
      //
      xrec = xlog;
      if(xtowc<0) xrec = 1-xlog;
      xrec -= 0.5;
      //
      myamp0 = 3.75-(3.75-2.75)/0.5*(v[2]+1);
      if(loge>5)
        myamp = myamp0;
      else
        myamp = myamp0*(1-(loge-5)*(loge-5)/10.);
      if(myamp<=0) continue;
      mysig = 0.09-(0.09-0.055)/0.5*(v[2]+1);
      if(mysig<=0) continue;
      par0 = 10.;
      par1 = myamp*exp(-v[i]*v[i]/2/mysig/mysig);
      //
      myxmax = fabs(1/par0*sqrt(par0/2/atan(par0/2)-1));
      mymax = fabs((myxmax-0.5/atan(par0*0.5)*atan(par0*myxmax)));
      delta = par1/mymax*(xrec-0.5/atan(par0*0.5)*atan(par0*xrec));
      //
      if(xtowc>0)
        pc[i] = p[i]-delta;
      else
        pc[i] = p[i]+delta;
    }

  Point corcentroid(pc[0],pc[1],pc[2]);

  return corcentroid;
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

