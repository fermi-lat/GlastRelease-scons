/** @mainpage package AcdDigi
* @author Heather Kelly
*
* @section intro Introduction
* This package contains all ACD digitization algorithms.
* The digitization takes as input the Monte Carlo data pertaining
* to the ACD volumes and determines the detector response of each
* PMT.
*
* Currently there are two algorithms for ACD digitization.
* AcdDigiAlg and AcdDigiOrgAlg.
*
* @section AcdDigiAlg
* Defines the current more realistic digitization routine for the ACD.
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
* @section AcdDigiOrgAlg
* Defines the first digitization routine for the ACD. 
* This version is left to allow easy comparison between this older version and the new version.
*
* <hr>
* @section jobOptions jobOptions
*
* <hr>
* @section notes release notes
* release.notes
* @section requirements requirements
* @include requirements
*
* @todo Fix up test program
*/

