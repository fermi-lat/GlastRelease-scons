// $Header$
// Mainpage for doxygen

/*! \mainpage package G4Generator

This package is the GLAST interface to Geant4. The steering is done with the main Algorithm, G4Generator.<br>

  For setup, it uses detModel to construct the detector, using a helper class 
  DetectorConstruction, which in turn uses G4SectionsVisitor and G4MaterialsVisitor.
  <br><br>
  Sensitive detector setup is now done with GlastDetectorManager, which creates the GlastDetector objects
  defined in the instrument class from XML files. It uses a vector of all the ids associated with 
  sensitive detectors to make the correspondence for passing hit information to the appropriate object.

  This will be gradually replaced or enhanced by saving the hit information itself on the TDS, for interpretation by 
  Algorithms in the GlastDigi package.
  <br><br>
  A special class DisplayManager uses the GuiSvc to display:
  - All sensitive detectors, ACD tiles, TKR planes, CAL logs , CAL diodes.
  - Ids of sensitive detectors
  - Hit detectors only
  - Steps
  - Tracks

  <br>
  The display of each can be enabled independently. 
  It will be relaced or supplemented by HepRep/WIRED.
  <br><br>
  Finally, the Algorithm IrfDisplay (copied from GismoGenerator) is used to display GlastDetector output, in the form
  of TdGlastData object on the TDS, as a test.

  A test program, under src/test, exercises everything. See the files below.

  <hr>
    \section requirements cmt/requirements
  \include requirements
    \section test src/test/jobOptions.txt
    \include src/test/jobOptions.txt
    \section test2 src/test/test_sources.xml
    \include src/test/test_sources.xml
*/

