/** @mainpage package AcdDigi
* @author Heather Kelly
*
* @section intro Introduction
* This package contains all ACD digitization algorithms.
* The digitization takes as input the Monte Carlo data pertaining
* to the ACD volumes and determines the detector response of each
* Photo Multiplier Tube (PMT).
*
* Currently there are three separate ACD digitization algorithms:
* AcdDigiAlg, AcdDigiMcIntHitAlg, and AcdDigiOrgAlg.
*
* @section AcdDigiAlg AcdDigiAlg
* Defines the current more realistic digitization routine for the ACD.
* Assumes ACD MC data is stored in McPositionHits
* - Optionally applies edge effects to adjust detected energy
* - Converts energy deposited in each ACD volume into Minimum Ionizing 
* Particles (MIPs)
* - MIPs are then converted to photoelectrons which are split evenly between 
* each PMT associated with an ACD volume.
* - Optionally adds Poisson Fluctuations to the number of photoelectrons seen 
* by each PMT.
* - Optionally adds Gaussian Noise for each readout:  PHA, veto, and CNO 
* discriminators.
* - Finds a PHA based on the number of MIPs associated with each PMT.
* - Determines how discriminator bits should be set on/off
*
* @section AcdDigiMcIntHitAlg AcdDigiMcIntHitAlg
* The first draft of a more realistic digitization routine for the ACD.
* Assumes ACD MC data is stored in McIntegratingHits.
* - Converts energy deposited in each ACD volume into MIPs
* - MIPs are then converted to photoelectrons which are split evenly between each PMT
* associated with an ACD volume.
* - Optionally adds Poisson Fluctuations to the number of photoelectrons seen 
* by each PMT.
* - Optionally adds Gaussian Noise for each readout:  PHA, veto and CNO 
* discriminators.
* - Finds a PHA based on the number of MIPs associated with each PMT.
* - Determines how discriminator bits should be set on/off
*
* @section AcdDigiOrgAlg AcdDigiOrgAlg
* Defines the first digitization routine for the ACD used for PDR studies. 
* This version is left to allow easy comparison between this older version and the 
* new versions.
*
* <hr>
* @section jobOptions jobOptions
* @param AcdDigiAlg.xmlFile
* The full path and filename of an input XML file containing constants for the 
* ACD digitization.
* @param AcdDigiAlg.autoCalibrate
* Boolean Flag that denotes whether or not to apply auto calibration
* 1 (One) denotes true, so that auto calibration will be applied
* 0 (zero) is false
* @param AcdDigiAlg.applyPoisson
* Boolean Flag denoting whether or not to apply Poisson fluctuations to the number of 
* photo electrons detected by a PMT.
* 1 (One) denotes true, so that Poisson fluctuations are applied.
* @param AcdDigiAlg.applyGaussianNoisse
* Boolean Flag denoting whether or not to apply gaussian noise to the electronics.
* Noise is applied to the PHA, veto and CNO disrciminators separately.
* 1 (One) denotes true, meaning that the noise is applied.
* @param AcdDigiAlg.edgeEffect
* Boolean Flag denoting whether or not to apply edge effects, where the position
* of a hit is used to determine how much energy was actually detected.
* 1 (One) denotes true, so that edge effects are taken into account.
* @param AcdDigiMcIntHitAlg.xmlFile
* The full path and filename of an input XML file containing constants for the
* ACD digitization.
* @param AcdDigiMcIntHitAlg.autoCalibrate
* Boolean Flag that denotes whether or not to apply auto calibration
* 1 (One) denotes true, so that auto calibration will be applied
* 0 (zero) is false
* @param AcdDigiMcIntHitAlg.applyPoisson
* Boolean Flag denoting whether or not to apply Poisson fluctuations to the number of 
* photo electrons detected by a PMT.
* 1 (One) denotes true, so that Poisson fluctuations are applied.
* @param AcdDigiMcIntHitAlg.applyGaussianNoisse
* Boolean Flag denoting whether or not to apply gaussian noise to the electronics.
* Noise is applied to the PHA, veto and CNO disrciminators separately.
* 1 (One) denotes true, meaning that the noise is applied.
* @param AcdDigiOrgAlg.xmlFile
* The full path and filename of an input XML file containing constants for 
* the ACD digitization.
*
* @section Tests Tests
* There are two tests associated with this package:  test_AcdDigi and test_AcdUtil.
* test_AcdUtil tests the AcdUtil class, making sure that its attempts to utilize
* the Poisson and Gaussian distributions in CLHEP are distributed as they should be.
*
* <hr>
* @section notes release notes
* release.notes
* @section requirements requirements
* @verbinclude requirements
*
* @todo Fix up test program
*/

