// $Header$
// Mainpage for doxygen

/** @mainpage package G4Generator
 *
 * @authors T.Burnett and R.Giannitrapani
 *
 * @section intro Introduction
 *
 * This package is the GLAST interface to Geant4. The steering is done with the
 * main Algorithm, G4Generator.<br>
 *
 * For setup, it uses GlastDetSvc to construct the detector, using a helper
 * class DetectorConstruction, which in turn uses the callback classes G4Media
 * and G4Geometry.  <br><br>
 *
 * Sensitive detector setup is now done with GlastDetectorManager. It uses a
 * vector of all the ids associated with sensitive detectors to make the
 * correspondence for passing hit information to the appropriate object.
 *
 * This will be gradually replaced or enhanced by saving the information itself
 * on the TDS, for interpretation by algorithms in the GlastDigi package.
 *
 * A special class DisplayManager uses the GuiSvc to display: 
 *          - All sensitive detectors, ACD tiles, TKR planes, CAL logs , CAL diodes.  
 *          - Ids of sensitive detectors 
 *          - Hit detectors only 
 *          - Steps 
 *          - Tracks
 *
 * <br> The display of each can be enabled independently.  It will be replaced or
 * supplemented by HepRep/WIRED.  <br>
 *
 * A test program, under src/test, exercises everything.
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr>
 * @section requirements requirements
 * @verbinclude requirements
 * <hr>
 * @todo Insert back the full detectors representation in the GUI from the
 * DetectorConstruction class; it needs some minor changes in detModel and
 * GlastDetSvc
 * @todo Check rotations sign in G4Geometry
 * @todo Fill McParticle tree
 *
 */

