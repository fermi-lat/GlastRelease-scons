// $Header$
// Mainpage for doxygen

/** @mainpage package CalRecon
 *
 * @authors A.Chekhtman, R.Terrier, J.A.Hernando 
 *
 * @section intro Introduction
 *
 *  CalRecon package reconstructs the energy and direction  of  incident particle from
 *   the calorimeter information. 
 *   
 *   Package contains 3 algorithms: CalXtalRecAlg, CalClustersAlg
 *   and CalDisplay.
 *
 *   CalXtalRecAlg takes the digitized calorimeter information from CalDigiCol
 *   as input, calculates the energy and position in each hitted crystal
 *   and stores this data into CalXtalRecCol.
 *    Energy for each crystal face is converted from adc value using simple
 *   linear formula:
 *   
 *    \f[
         E = E_{max} * (ADC - PED)/(ADC_{max} - PED)
     \f]   
 *   where \f$ E_{max} \f$ - maximum energy for used energy range,
 *  \f$ PED \f$ - pedestal, \f$ ADC_{max} \f$ - maximum ADC value.
 *   These constants are the same for all crystals and are defined
 *   in xml file flightCALResponse.xml .
 *
 *   Position along the crystal direction is calculated from asymmetry between
 *   energies reconstructed from positive and negative faces of the crystal
 *
 *   \f[ pos = \frac{E_{pos} - E_{neg}}{E_{pos} + E_{neg}} * 
               \frac{1+lightAtt}{1-lightAtt} * \frac{L_{crystal}}{2}
     \f]
 *   where \f$ L_{crystal} \f$ - crystal length and \f$ lightAtt \f$ - 
 *   the ratio of signal from "far" crystal face to the signal
 *   from"near" crystal face in the case when the energy deposition is close
 *   to one end of the crystal.
 *  
 *   CalClustersAlg calculates the energy, position and direction for
 *   calorimeter clusters and applies energy corrections.
 *   Actually there is no real clustering algorithm implemented,
 *   CalClusterAlg now considers, that there is only one cluster including
 *   all hitted calorimeter crystals. For this "cluster" the algorithm
 *   calculates the following parameters and stores them in the 
 *   Event::CalClusterCol object:
 *     - total energy sum 
 *     - energy corrected by profile fitting method
 *     - energy corrected by last layer correlation method
 *     - energy sum for each layer
 *     - average position 
 *     - average position for each layer
 *     - direction of incident particle
 *    
 *
 *   CalDisplay algorithm provides the display of reconstructed data.
 *
 *  
 * @todo modify CalXtalRecAlg to use real calibration data
 *
 * @todo implement real clustering in CalClustersAlg to determine the 
 *       energies of electron and positron produced
 *       by low energy photon
 *
 * @todo move out hardwired constants describing the energy correction fit
 *       from CalClustersAlg::Leak function into a data file;
 *       update this constants to reflect the actual detector geometry
 *        
 * @todo implement the low energy correction based on number of tracker
 *       hits 
 *      
 * <hr>
 * @section notes release notes
 * @include release.notes
 * @section requirements requirements
 * @verbinclude requirements

*/

