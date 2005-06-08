
// File and Version information:
// $Header$
//
//  Implementation file of Rot clone for CalParams
//  
// Authors:
//
//    Tracy Usher
//    David Chamont
//

#include <reconRootData/CalParams.h>
#include <iostream>

ClassImp(CalParams)

CalParams::CalParams
( double energy, double eneError,
  TVector3 centroid,
  double cntdxx, double cntdxy, double cntdxz, double cntdyy, double cntdyz, double cntdzz,
  TVector3 axis,
  double axsdxx, double axsdxy, double axsdxz, double axsdyy, double axsdyz, double axsdzz )
{
    init
    ( energy,eneError,
      centroid,cntdxx,cntdxy,cntdxz,cntdyy,cntdyz,cntdzz,
      axis,axsdxx,axsdxy,axsdxz,axsdyy,axsdyz,axsdzz ) ;
}

void CalParams::Clear( Option_t * )
{
    init
    ( 0.,0.,
      TVector3(0.,0.,0.),0.,0.,0.,0.,0.,0.,
      TVector3(0.,0.,0.),0.,0.,0.,0.,0.,0. ) ;
}

void CalParams::init
( Double_t energy, Double_t eneError,
  TVector3 centroid,
  Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
  Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
  TVector3 axis,
  Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
  Double_t axsdyy, Double_t axsdyz, Double_t axsdzz )
{
    m_energy = energy ;
    m_eneError = eneError ;
    m_clusterCentroid = centroid ;
    m_cenxx = cntdxx ;
    m_cenxy = cntdxy ;
    m_cenxz = cntdxz ;
    m_cenyy = cntdyy ;
    m_cenyz = cntdyz ;
    m_cenzz = cntdzz ;
    m_clusterAxis = axis ;
    m_axisxx = axsdxx ;
    m_axisxy = axsdxy ;
    m_axisxz = axsdxz ;
    m_axisyy = axsdyy ;
    m_axisyz = axsdyz ;
    m_axiszz = axsdzz ;
}

void CalParams::Print( Option_t * ) const
{
    std::cout
      << m_energy << " " << m_eneError << "\n"
      << m_clusterCentroid.X() << " " << m_clusterCentroid.Y() << " " << m_clusterCentroid.Z() << "\n"
      << m_cenxx  << " " << m_cenxy << " " << m_cenxz << "\n"
      << m_cenxy  << " " << m_cenyy << " " << m_cenyz << "\n"
      << m_cenxz  << " " << m_cenyz << " " << m_cenzz << "\n"
      << m_clusterAxis.X() << " " << m_clusterAxis.Y() << " " << m_clusterAxis.Z() << "\n"
      << m_axisxx  << " " << m_axisxy << " " << m_axisxz << "\n"
      << m_axisxy  << " " << m_axisyy << " " << m_axisyz << "\n"
      << m_axisxz  << " " << m_axisyz << " " << m_axiszz ;
}

