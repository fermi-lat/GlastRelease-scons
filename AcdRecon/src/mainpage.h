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
 * contains the raw ACD event data as well as the geometrical representation of the ACD.
 *
 * TDS Inputs:
 *  - Event::AcdDigiCol: all the AcdDigi objects
 *  - Event::McParticleCol: MC tracks, obviously for MC only
 *
 * TDS Outputs:
 *  - Event::AcdRecon object, which contains
 *    - Event::AcdTkrIntersectionCol: all the intersections with GEANT model 
 *    - Event::AcdTkrHitPocaCol all the POCA calculations
 *    - Event::AcdTkrGapPocaCol: POCA w.r.t. gaps and ribbons in the GEANT model 
 *    - Event::AcdTkrPointCol: Data about where the tracks exits the ACD volume
 *    - Event::AcdSplashVarsCol: NOT Filled!! (Data about backsplash from tracks)
 *    - Some numbers we extracted during the recon process
 *      - Number of hit tiles
 *      - Number of hit ribbons
 *      - Max active distance over all tiles, ID of tile w/ max active distance
 *      - Max active distance over all ribbons, ID of ribbon w/ max active distance
 *      - Max active distance for all rows of tiles of ACD
 *      - Minimum distance to a corner ray
 *    - Some MC quantaties
 *      - Total tile energy (MC)
 *      - Total ribbon energy (MC)
 *      - twin vectors of AcdID and MC energy
 *    - Some useless crap
 *      - Gamma DOCA (Deprecated!!)
 *      - Min DOCA over all tracks, ID of tile w/ min DOCA (Deprecated!!)
 *
 * @section jobOptions jobOptions
 * 
 * @subsection AcdPha2MipTool_JO AcdPha2MipTool
 *  - AcdCalibSvc ["AcdCalibSvc"]  : Name of Acd Calibration SVC to use
 *  - PHATileCut [0.]              : Ignore all tiles with pedestal subtracted PHA below this value
 *  - MIPSTileCut [0.]             : Ignore all tiles with MIP equivalent below this value
 *  - PHARibbonCut [0.]            : Ignore all ribbons with pedestal subtracted PHA below this valu
 *  - MIPSRibbonCut [0.]           : Ignore all ribbons with MIP equivalent below this value 
 *
 * @subsection AcdPocaTool_JO AcdPocaTool
 *  - distanceCut [1999.]   : Filter out POCA when doca is > this value
 *  - sigmaCut [5.]         : Unused!!  (Filter out POCA when doca/docaError) is > this value
 *
 * @subsection AcdReconAlg _JO AcdReconAlg
 *  - Tool names:
 *    -  intersectionToolName["AcdTkrIntersectTool"]
 *    -  hitToolName["AcdPha2MipTool"]
 *    -  pocaToolName["AcdPocaTool"]
 *    -  propToolName["G4PropagationTool"]
 *  - doBackSplash[false] : do turn on backsplash calculations (caveat emptor)
 *
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

