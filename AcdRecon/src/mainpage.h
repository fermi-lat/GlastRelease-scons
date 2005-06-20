/** @mainpage package AcdRecon
 * @author Heather Kelly
 *
 * @section intro Introduction
 * This package contains the ACD (Anti-Coincidence Detector) reconstruction
 * algorithms.  The primary class in this package is AcdReconAlg which is a 
 * Gaudi Algorithm.
 *
 * @section AcdReconAlg AcdReconAlg
 * A Gaudi algorithm that uses the Event::AcdDigi collection on the TDS, which
 * contains the ACD detector information.  AcdReconAlg is based on the original VetoRecon 
 * algorithm available within glastsim.
 *
 * The primary output from the ACD reconstruction is the EventModel::AcdRecon object
 * stored on the TDS.  AcdReconAlg computes the following:
 * - Total energy deposited in the ACD detectors (counting only those above 
 * veto threshold)
 * - Total number of ACD tiles above veto threshold
 * - Minimum Distance of Closest Approach (DOCA) 
 * - A list of DOCA values containing the minimum DOCA for the top and side rows
 * - Minimum Active Distance
 * - A list of Active Distance values containing the minimum values for the 
 * top and side rows
 * The DOCA and Active Distance quantities are computed using the ACD detector hits
 * and the TkrRecon reconstructed track collection.  DOCA is calculated by finding
 * the minimum distance between the center of hit ACD tiles and all found tracks.
 * Active Distance is calculated by finding the minimum distance between the edge
 * of hit ACD tiles and all found tracks.
 * 
 * @section jobOptions jobOptions
 * No jobOptions are currently utilized in the AcdRecon package.
 * 
 <hr>
 * @section jobOptions jobOptions
 * NONE
 *
 * @section Tests Tests
 * There is one test routine available:  test_AcdRecon.
 * This test uses a digi.root file that stores the AcdDigi values to be used when running
 * AcdReconAlg.  Hence, RootIo is required to run this test.
 *
 * @section notes release.notes
 * release.notes
 * <hr> 
 * @section requirements requirements
 * @verbinclude requirements
 * <hr>
 */

