/** @mainpage package AcdRecon
 * @author Heather Kelly
 *
 * @section intro Introduction
 * This package contains the ACD (Anti-Coincidence Detector) reconstruction
 * algorithms.
 * The algorithms are implemented as Gaudi algorithms.
 * The primary output from the ACD reconstruction is the EventModel::AcdRecon object
 * stored on the TDS.
 * AcdReconAlg computes the following:
 * - Total number of ACD tiles above veto threshold
 * - Minimum Distance of Closest Approach
 * - A list of DOCA values for the top and side rows
 * - Minimum Active Distance
 * - A list of Active Distance values for the top and side rows
 * 
 * @section requirements requirements
 * @include requirements
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr> 
 */

