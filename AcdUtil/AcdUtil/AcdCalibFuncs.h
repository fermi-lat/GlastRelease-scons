#ifndef ACDCALIBFUNCS_H
#define ACDCALIBFUNCS_H


#include "GaudiKernel/StatusCode.h"

namespace AcdCalib {


  // Estimate the number of photo electrons seen in the PMT
  // This is the nominal number, ie, the mean of a Possion distribution
  StatusCode photoElectronsFromEnergy(const double& E, const double& pe_per_MeV, double& pe);

  // Convert the number of observed photo electrons into and signal
  // The signal is express in mip equivalent units
  StatusCode lightYeildMipEquivalent(const double& pe, const double& pe_per_mip, double& signal_mips);
  
  // Determine the low range PHA value corrsponding to a particular MIP equivalent
  StatusCode PHA_lowRange(const double& signalMips, const double& pedestal, const double& mipPeak, 
			  unsigned short& PHA);

  // Determine the high range PHA value corrsponding to a particular MIP equivalent
  StatusCode PHA_highRange(const double& signalMips, const double& pedestal, const double& slope, const double& saturation, 
			   unsigned short& PHA);
  
  // Determine the MIP equivalent value of a particular low range PHA value
  StatusCode mipEquivalent_lowRange(unsigned short PHA, const double& pedestal, const double& mipPeak, double& mips);

  // Determine the MIP equivalent value of a particular high range PHA value
  StatusCode mipEquivalent_highRange(unsigned short PHA, const double& pedestal, const double& slope, const double& saturation, double& mips);

  // Determine the Z value of a particular charge deposit
  StatusCode Z_value(const double& mips, const double& pathFactor, double& Z);
 
}


#endif
