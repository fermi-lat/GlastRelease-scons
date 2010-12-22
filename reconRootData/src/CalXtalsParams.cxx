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
    "Raw sum of xtal energies = " << m_xtalRawEneSum << "\n" <<
    "Corrected sum of xtal energies = " << m_xtalCorrEneSum << "\n" <<
    "Rms of xtal energy distribution = " << m_xtalEneRms << "\n" <<
    "Skewness of xtal energy distribution = " << m_xtalEneSkewness << std::endl;
}

// dummy data, just for tests
void CalXtalsParams::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ )
 {
   init(1,2,3,4.,5.,6.,7.);
 }

Bool_t CalXtalsParams::CompareInRange( const CalXtalsParams & cp, const std::string & name ) const
{  
  Bool_t result = true ;
  // Not implemented, yet.
  std::cout << "CalXtalsParams::CompareInRange() not implemented, yet." << std::endl;
  return result;
}

