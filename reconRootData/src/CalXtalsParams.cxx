//  Implementation file of Root clone for CalXtalsParams
//  
// Authors:
//
//    Luca Baldini
//    Johan Bregeon
//

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
  Double_t xtalEneRms,
  Double_t xtalEneSkewness )
{
  init
    ( numXtals,numTruncXtals,numSaturatedXtals,
      xtalRawEneSum,xtalCorrEneSum,xtalEneRms,xtalEneSkewness );
}

void CalXtalsParams::init(Int_t numXtals, Int_t numTruncXtals, Int_t numSaturatedXtals,
			  Double_t xtalRawEneSum, Double_t xtalCorrEneSum,
			  Double_t xtalEneRms, Double_t xtalEneSkewness)
{
  m_numXtals          = numXtals;
  m_numTruncXtals     = numTruncXtals;
  m_numSaturatedXtals = numSaturatedXtals;
  m_xtalRawEneSum     = xtalRawEneSum;
  m_xtalCorrEneSum    = xtalCorrEneSum;
  m_xtalEneRms        = xtalEneRms;
  m_xtalEneSkewness   = xtalEneSkewness;
}

void CalXtalsParams::Clear( Option_t * )
{
  m_numXtals          = 0;
  m_numTruncXtals     = 0;
  m_numSaturatedXtals = 0;
  m_xtalRawEneSum     = 0.;
  m_xtalCorrEneSum    = 0.;
  m_xtalEneRms        = -1.;
  // Find a sensible value, here.
  m_xtalEneSkewness   = 0.;
}

void CalXtalsParams::Print( Option_t * ) const
{
  std::cout <<
    "Number of xtals = " << m_numXtals << "\n" <<
    "Truncated number of xtals = " << m_numTruncXtals << "\n" <<
    "Number of saturated xtals = " << m_numSaturatedXtals << "\n" <<
    "Raw sum of xtal energies = " << m_xtalRawEneSum << " MeV\n" <<
    "Corrected sum of xtal energies = " << m_xtalCorrEneSum << " MeV\n" <<
    "Rms of xtal energy distribution = " << m_xtalEneRms << " MeV\n" <<
    "Skewness of xtal energy distribution = " << m_xtalEneSkewness << std::endl;
}

// dummy data, just for tests
void CalXtalsParams::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ )
 {
   init(1,2,3,4.,5.,6.,7.);
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
  result = rootdatautil::CompareInRange(getXtalEneRms(),xp.getXtalEneRms(),
					"Rms of xtal energy distribution") && result ;
  result = rootdatautil::CompareInRange(getXtalEneSkewness(),xp.getXtalEneSkewness(),
					"Skewness of xtal energy distribution") && result ;
  
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
