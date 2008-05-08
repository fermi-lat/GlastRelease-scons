#ifndef ACDCALIBFUNCS_H
#define ACDCALIBFUNCS_H

// Gaudi, facilities, interfaces
#include "GaudiKernel/StatusCode.h"

/**
 * @brief Some functions that do conversion using ACD calibrations
 *
 **/

namespace AcdCalib {


  /**
   * @brief Estimate the number of photo-electrons seen in the PMT
   *
   * This is the nominal number, ie, the mean of a Possion distribution
   *
   * @param E energy deposited according to GEANT (and correct for edge effects)
   * @param pe_per_MeV conversion factor from energy to p.e.
   * @param pe mean number of photo-electrons seen at this energy
   * @return StatusCode indicated success or failure of conversion
   */
  StatusCode photoElectronsFromEnergy(const double& E, const double& pe_per_MeV, double& pe);

  /**
   * @brief Convert the number of observed photo-electrons into a MIP-equivalent signal
   *
   * The signal is express in mip equivalent units.  
   *
   * @param pe mean number of photo-electrons seen at this energy
   * @param pe_per_mip conversion factor photo-electrons to MIPS
   * @param signal_mips signal size in mip equivalent units
   * @return StatusCode indicated success or failure of conversion
   */
  StatusCode lightYeildMipEquivalent(const double& pe, const double& pe_per_mip, double& signal_mips);
  
  /**
   * @brief Determine the low range PHA value corrsponding to a particular MIP equivalent
   *
   * This just uses a linear transformation:  
   *     PHA = pedestal + ( signalMips * mipPeak )
   *
   * @param signal_mips signal size in mip equivalent units
   * @param pedestal in PHA units, 
   * @param mipPeak in PHA units above pedestal
   * @param PHA signal expressed in PHA 
   * @return StatusCode indicated success or failure of conversion
   */
  StatusCode PHA_lowRange(const double& signalMips, const double& pedestal, const double& mipPeak, 
			  unsigned short& PHA);

  /**
   * @brief Determine the high range PHA value corrsponding to a particular MIP equivalent
   *
   * This uses a transformation that goes from linear to saturating:
   *     PHA = pedestal + ( ( signalMips * slope * saturation ) / ( ( signalMips * slope ) + ( saturation ) ) )
   *
   * @param signal_mips signal size in mip equivalent units
   * @param pedestal in PHA units, 
   * @param slope conversion for mip to PHA for low values
   * @param saturation saturations for in PHA 
   * @param PHA signal expressed in PHA 
   * @return StatusCode indicated success or failure of conversion
   */
  StatusCode PHA_highRange(const double& signalMips, const double& pedestal, const double& slope, const double& saturation, 
			   unsigned short& PHA);
  
  /**
   * @brief Determine the MIP equivalent value of a particular low range PHA value
   *
   * This just uses a linear transformation:  
   *     mips = (PHA - pedestal) / mipPeak
   *
   * @param PHA signal expressed in PHA 
   * @param pedestal in PHA units, 
   * @param mipPeak in PHA units above pedestal
   * @param mips signal size in mip equivalent units
   * @return StatusCode indicated success or failure of conversion
   */
  StatusCode mipEquivalent_lowRange(unsigned short PHA, const double& pedestal, const double& mipPeak, double& mips);

  /**
   * @brief Determine the MIP equivalent value of a particular high range PHA value
   *
   * This uses a transformation that goes from linear to saturating:
   *     mips = (PHA - pedestal) * slope * saturation / ( saturation + ( slope * (PHA - pedestal) ) )
   *
   * @param PHA signal expressed in PHA 
   * @param pedestal in PHA units, 
   * @param slope conversion for mip to PHA for low values
   * @param saturation saturations for in PHA 
   * @param mips signal size in mip equivalent units
   * @return StatusCode indicated success or failure of conversion
   */
  StatusCode mipEquivalent_highRange(unsigned short PHA, const double& pedestal, const double& slope, const double& saturation, double& mips);
  /**
   * @brief Determine the change in pedestal b/c of the coherent noise
   *
   * This describes ringing, ie, an falling exponential X an oscillation
   *     deltaPed = amplitude * exp( - gemDT * decay ) * sin ( ( gemDT * freq ) + phase )
   *
   * @param gemDT is the "GemDeltaEventTime", ie time since last event
   * @param amplitude of the pedestal oscillations 
   * @param decay of the exponential fall off
   * @param freq. of the oscillations
   * @param phase of the oscillations
   * @param deltaPed change in pedestal 
   * @return StatusCode indicated success or failure of conversion
   */
  StatusCode coherentNoise(unsigned gemDT, 
			   const double& amplitude, const double& decay, const double& freq, const double& phase,  
			   double& deltaPed );

  /**
   * @brief Determine the Z value of a particular charge deposit
   *
   * FIMXE.  This is broken
   *
   * @param mips signal size in mip equivalent units
   * @param pathFactor is ratio of pathlength to tile thickness 
   * @param Z is the Z value (estimated atomic number)
   * @return StatusCode indicated success or failure of conversion
   */
  StatusCode Z_value(const double& mips, const double& pathFactor, double& Z);
 
}


#endif
