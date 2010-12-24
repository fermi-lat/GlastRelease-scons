
// File and Version information:
// $Header$
//
//  Implementation file of Rot clone for CalFitParams
//  
// Authors:
//
//    Tracy Usher
//    David Chamont
//    Luca Baldini
//

#include <reconRootData/CalFitParams.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalFitParams)

CalFitParams::CalFitParams
( Double_t energy, Double_t eneError,
  TVector3 centroid,
  Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
  Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
  TVector3 axis,
  Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
  Double_t axsdyy, Double_t axsdyz, Double_t axsdzz,
  Int_t numFitLayers, Double_t chiSquare )
{
  init
    ( energy,eneError,
      centroid,cntdxx,cntdxy,cntdxz,cntdyy,cntdyz,cntdzz,
      axis,axsdxx,axsdxy,axsdxz,axsdyy,axsdyz,axsdzz,
      numFitLayers, chiSquare ) ;
}

CalFitParams::CalFitParams
( int numFitLayers, double chiSquare,
  TVector3 centroid,
  double cntdxx, double cntdxy, double cntdxz,
  double cntdyy, double cntdyz, double cntdzz,
  TVector3 axis,
  double axsdxx, double axsdxy, double axsdxz,
  double axsdyy, double axsdyz, double axsdzz )
{
  init
    ( numFitLayers, chiSquare,
      centroid,cntdxx,cntdxy,cntdxz,cntdyy,cntdyz,cntdzz,
      axis,axsdxx,axsdxy,axsdxz,axsdyy,axsdyz,axsdzz ) ;
}

void CalFitParams::init
( Double_t energy, Double_t eneError,
  TVector3 centroid,
  Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
  Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
  TVector3 axis,
  Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
  Double_t axsdyy, Double_t axsdyz, Double_t axsdzz,
  Int_t numFitLayers, Double_t chiSquare )
{
  setEnergy(energy);
  setEnergyErr(eneError);
  setCentroid(centroid);
  setxPosxPos(cntdxx);
  setxPosyPos(cntdxy);
  setxPoszPos(cntdxz);
  setyPosyPos(cntdyy);
  setyPoszPos(cntdyz);
  setzPoszPos(cntdzz);
  setAxis(axis);
  setxDirxDir(axsdxx);
  setxDiryDir(axsdxy);
  setxDirzDir(axsdxz);
  setyDiryDir(axsdyy);
  setyDirzDir(axsdyz);
  setzDirzDir(axsdzz);
  m_nFitLayers = numFitLayers ;
  m_chiSquare  = chiSquare ;
}

void CalFitParams::init
( Int_t numFitLayers, Double_t chiSquare,
  TVector3 centroid,
  Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
  Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
  TVector3 axis,
  Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
  Double_t axsdyy, Double_t axsdyz, Double_t axsdzz )
{
  setEnergy(-1.);
  setEnergyErr(-1.);
  setCentroid(centroid);
  setxPosxPos(cntdxx);
  setxPosyPos(cntdxy);
  setxPoszPos(cntdxz);
  setyPosyPos(cntdyy);
  setyPoszPos(cntdyz);
  setzPoszPos(cntdzz);
  setAxis(axis);
  setxDirxDir(axsdxx);
  setxDiryDir(axsdxy);
  setxDirzDir(axsdxz);
  setyDiryDir(axsdyy);
  setyDirzDir(axsdyz);
  setzDirzDir(axsdzz);
  m_nFitLayers = numFitLayers ;
  m_chiSquare  = chiSquare ;
}

void CalFitParams::Clear( Option_t * )
{
  init
    ( 0.,0.,
      TVector3(0.,0.,0.),0.,0.,0.,0.,0.,0.,
      TVector3(0.,0.,0.),0.,0.,0.,0.,0.,0.,
      0,0. ) ;
}

void CalFitParams::Print( Option_t * ) const
{
  // Call the base class method...
  CalParams::Print();
  // ... then print the additional stuff.
  std::cout << "\n" <<
    "Number of layers for the fit: " << m_nFitLayers << "\n" <<
    "Fit chisquare: " << m_chiSquare << std::endl;
}

// dummy data, just for tests
void CalFitParams::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ )
 {
    init
    ( 1,2.,
      TVector3(3.,4.,5.),6.,7.,8.,9.,10.,11.,
      TVector3(12.,13.,14.),15.,16.,17.,18.,19.,20.) ;
 }

Bool_t CalFitParams::CompareInRange( const CalFitParams & cp, const std::string & name ) const {
    
    Bool_t result = true ;
    
    result = rootdatautil::CompareInRange(getFitLayers(),cp.getFitLayers(),"Fit Layers") && result ;
    result = rootdatautil::CompareInRange(getChiSquare(),cp.getChiSquare(),"ChiSquare") && result ;
    result = rootdatautil::CompareInRange(getCentroid(),cp.getCentroid(),"Centroid") && result ;
    result = rootdatautil::CompareInRange(getxPosxPos(),cp.getxPosxPos(),"xPosxPos") && result ;
    result = rootdatautil::CompareInRange(getxPosyPos(),cp.getxPosyPos(),"xPosyPos") && result ;
    result = rootdatautil::CompareInRange(getxPoszPos(),cp.getxPoszPos(),"xPoszPos") && result ;
    result = rootdatautil::CompareInRange(getyPosyPos(),cp.getyPosyPos(),"yPosyPos") && result ;
    result = rootdatautil::CompareInRange(getyPoszPos(),cp.getyPoszPos(),"yPoszPos") && result ;
    result = rootdatautil::CompareInRange(getzPoszPos(),cp.getzPoszPos(),"zPoszPos") && result ;
    result = rootdatautil::CompareInRange(getAxis(),cp.getAxis(),"Axis") && result ;
    result = rootdatautil::CompareInRange(getxDirxDir(),cp.getxDirxDir(),"xDirxDir") && result ;
    result = rootdatautil::CompareInRange(getxDiryDir(),cp.getxDiryDir(),"xDiryDir") && result ;
    result = rootdatautil::CompareInRange(getxDirzDir(),cp.getxDirzDir(),"xDirzDir") && result ;
    result = rootdatautil::CompareInRange(getyDiryDir(),cp.getyDiryDir(),"yDiryDir") && result ;
    result = rootdatautil::CompareInRange(getyDirzDir(),cp.getyDirzDir(),"yDirzDir") && result ;
    result = rootdatautil::CompareInRange(getzDirzDir(),cp.getzDirzDir(),"zDirzDir") && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for CalFitParams"<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;
        
}

