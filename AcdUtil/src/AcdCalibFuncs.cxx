#include "../AcdUtil/AcdTileFuncs.h"

#include "CLHEP/Matrix/Matrix.h"

namespace AcdCalib {
 
  // Estimate the number of photo electrons seen in the PMT
  // This is the nominal number, ie, the mean of a Possion distribution
  StatusCode photoElectronsFromEnergy(const double& E, const double& pe_per_MeV, double& pe) {
    pe = E * pe_per_MeV;
    return StatusCode::SUCCESS;
  }

  // Convert the number of observed photo electrons into and signal
  // The signal is express in mip equivalent units
  StatusCode lightYeildMipEquivalent(const double& pe, const double& pe_per_mip, double& signal_mips) {
    if ( pe_per_mip <= 0. ) {
      signal_mips = 0.;
      return StatusCode::FAILURE;
    }
    signal_mips = pe / pe_per_mip;
    return StatusCode::SUCCESS;
  }
  
  // Determine the low range PHA value corrsponding to a particular MIP equivalent
  StatusCode PHA_lowRange(const double& signalMips, const double& pedestal, const double& mipPeak, unsigned short& PHA) {
    double PHA_f = signalMips * mipPeak;
    PHA_f += pedestal;
    PHA = PHA_f < 0 ? 0 : (unsigned short)PHA_f;
    return StatusCode::SUCCESS;
  }

  // Determine the high range PHA value corrsponding to a particular MIP equivalent
  StatusCode PHA_highRange(const double& signalMips, const double& pedestal, const double& slope, const double& saturation, unsigned short& PHA) {
    double top = signalMips * slope * saturation;
    double bot = saturation + ( signalMips * slope );
    if ( bot <= 0. ) {
      PHA = 0;
      return StatusCode::FAILURE;
    }
    double PHA_f = top/bot;
    PHA_f += pedestal;
    PHA = PHA_f < 0 ? 0 : (unsigned short)PHA_f;
    return StatusCode::SUCCESS;
  }
  
  // Determine the MIP equivaltent value of a particular low range PHA value
  StatusCode mipEquivalent_lowRange(unsigned short  PHA, const double& pedestal, const double& mipPeak, double& mips) {
    if ( mipPeak <= 0. ) {
      mips = 0.;
      return StatusCode::FAILURE;
    }
    double PHApedReduced = (double)PHA - pedestal;
    mips = PHApedReduced / mipPeak;
    return StatusCode::SUCCESS;
  }

  // Determine the MIP equivaltent value of a particular high range PHA value
  StatusCode mipEquivalent_highRange(unsigned short PHA, const double& pedestal, const double& slope, const double& saturation, double& mips) {
    double PHApedReduced = (double)PHA - pedestal;
    double top = PHApedReduced * saturation;
    double bot = ( saturation - PHApedReduced ) * slope;
    if ( bot <= 0 ) {
      mips = 0;
      return StatusCode::FAILURE;
    }
    mips = top/bot;
    return StatusCode::SUCCESS;
  }

  // Determine the Z value of a particular charge deposit
  StatusCode Z_value(const double& mips, const double& pathFactor, double& Z) {
    if ( pathFactor <= 0. ) {
      Z = 0;
      return StatusCode::FAILURE;
    }
    double mipsReduced = mips/pathFactor;
    return StatusCode::SUCCESS;
  }

}
