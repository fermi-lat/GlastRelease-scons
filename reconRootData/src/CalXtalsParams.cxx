
/**
   @file CalXtalsParams.cxx

   @brief Implementation of the CalXtalsParams class.
      
   @author Luca Baldini (luca.baldini@pi.infn.it)
   
   $Revision$
   $Date$
   $Header$
*/

#include <reconRootData/CalXtalsParams.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalXtalsParams)

CalXtalsParams::CalXtalsParams()
{
  Clear();
}

CalXtalsParams::CalXtalsParams
( Int_t numXtals,
  Int_t numTruncXtals,
  Int_t numSaturatedXtals,
  Double_t xtalRawEneSum,
  Double_t xtalCorrEneSum,
  Double_t xtalEneMax,
  Double_t xtalEneRms,
  Double_t xtalEneSkewness,
  TVector3 centroid )
{
  init
    ( numXtals,numTruncXtals,numSaturatedXtals,xtalEneMax,
      xtalRawEneSum,xtalCorrEneSum,xtalEneRms,xtalEneSkewness,centroid );
}

void CalXtalsParams::Clear( Option_t * )
{
  init
    ( 0.,0.,0.,0.,0.,0.,-1.,0.,TVector3(0.,0.,0.) ) ;
}

void CalXtalsParams::init(Int_t numXtals, Int_t numTruncXtals, Int_t numSaturatedXtals,
			  Double_t xtalRawEneSum, Double_t xtalCorrEneSum, Double_t xtalEneMax,
			  Double_t xtalEneRms, Double_t xtalEneSkewness, TVector3 centroid)
{
  m_numXtals          = numXtals;
  m_numTruncXtals     = numTruncXtals;
  m_numSaturatedXtals = numSaturatedXtals;
  m_xtalRawEneSum     = xtalRawEneSum;
  m_xtalCorrEneSum    = xtalCorrEneSum;
  m_xtalEneMax        = xtalEneMax;
  m_xtalEneRms        = xtalEneRms;
  m_xtalEneSkewness   = xtalEneSkewness;
  m_centroid          = centroid;
}

void CalXtalsParams::Print( Option_t * ) const
{
  std::cout <<
    "Number of xtals = " << m_numXtals << "\n" <<
    "Truncated number of xtals = " << m_numTruncXtals << "\n" <<
    "Number of saturated xtals = " << m_numSaturatedXtals << "\n" <<
    "Raw sum of xtal energies = " << m_xtalRawEneSum << " MeV\n" <<
    "Corrected sum of xtal energies = " << m_xtalCorrEneSum << " MeV\n" <<
    "Maximum xtal energy = " << m_xtalEneMax << " MeV\n" <<
    "Rms of xtal energy distribution = " << m_xtalEneRms << " MeV\n" <<
    "Skewness of xtal energy distribution = " << m_xtalEneSkewness <<
    "Centroid = (" << m_centroid.x() << ", " << m_centroid.y() << ", "
			 << m_centroid.z() << ") mm" << std::endl;
}

// dummy data, just for tests
void CalXtalsParams::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ )
 {
   init(1,2,3,4.,5.,5.5,6.,7.,TVector3(-1.,-1.,-1.));
 }

Bool_t CalXtalsParams::CompareInRange( const CalXtalsParams & xp, const std::string & name ) const
{  
  Bool_t result = true ;
  result = rootdatautil::CompareInRange(getNumXtals(),xp.getNumXtals(),
					"Number of xtals") && result ;
  result = rootdatautil::CompareInRange(getNumTruncXtals(),xp.getNumTruncXtals(),
					"Truncated number of xtals") && result ;
  result = rootdatautil::CompareInRange(getNumSaturatedXtals(),xp.getNumSaturatedXtals(),
					"Number of saturated xtals") && result ;
  result = rootdatautil::CompareInRange(getXtalRawEneSum(),xp.getXtalRawEneSum(),
					"Raw xtal energy sum") && result ;
  result = rootdatautil::CompareInRange(getXtalCorrEneSum(),xp.getXtalCorrEneSum(),
					"Corrected xtal energy sum") && result ;
  result = rootdatautil::CompareInRange(getXtalEneMax(),xp.getXtalEneMax(),
					"Maximum xtal energy") && result ;
  result = rootdatautil::CompareInRange(getXtalEneRms(),xp.getXtalEneRms(),
					"Rms of xtal energy distribution") && result ;
  result = rootdatautil::CompareInRange(getXtalEneSkewness(),xp.getXtalEneSkewness(),
					"Skewness of xtal energy distribution") && result ;
  result = rootdatautil::CompareInRange(getCentroid(),xp.getCentroid(),
					"Centroid of xtal collection") && result ;
  
  if (!result) {
    if ( name == "" ) {
      std::cout<<"Comparison ERROR for CalXtalsParams"<<std::endl ;
    }
    else {
      std::cout<<"Comparison ERROR for "<<name<<std::endl ;
    }
  }
  
  return result;
}
