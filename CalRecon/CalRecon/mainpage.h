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
 *   
 *   CalClustersAlg calculates the energy position and direction for calorimeter clusters
 *   and applies energy corrections.
 *    
 *   CalDisplay algorithm provides the display of reconstructed data.
 *
 *  
 *    
 *      
 * <hr>
 * @section notes release notes
 * @include release.notes
 * @section requirements requirements
 * @verbinclude requirements

*/

