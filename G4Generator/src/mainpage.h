// $Header$
// Mainpage for doxygen

/** @mainpage package G4Generator
 *
 * @authors T.Burnett and R.Giannitrapani
 *
 * @section description Description
 *
 * This package is the GLAST interface to the Montecarlo simulations toolkit
 * Geant4 (see <a href="http://geant4.web.cern.ch">here</a> 
 * for more information on this C++ toolkit). 
 * The steering is done with the main Algorithm, G4Generator.<br>
 * 
 * The main simulation is controlled by RunManager, a customized and stripped
 * version of the Geant4 standard G4RunManager singleton. Since GLAST main event
 * loop is driven by GAUDI and it will not use any graphics or data persistency
 * feature of Geant4, we have included in the RunManager only the real necessary
 * part for setup and run the generator. RunManager itself uses the following
 * classes 
 *
 * - DetectorConstruction: this class, derived by a standard
 * G4VUserDetectorConstruction, setup both the material list and the geometry of
 * the detector. In our case these information are stored in xml files; to
 * access them the DetectorConstruction class uses methods of the Gaudi service
 * GlastDetSvc and two classes, G4Geometry and G4Media, that implements abstract
 * interfaces of the GlastSvc package. All this is based on a Visitor mechanism
 * that separates the concrete implementation of the xml geometry files to any
 * possible client (simulation, reconstruction, event display etc.). G4Geometry
 * implements methods to traverse the geometry of GLAST and build a concrete
 * Geant4 representation of it; G4Media do the same for the materials
 * definitions.
 * - PhysicsList: this class, derived by a standard G4VModularPhysicsList, is
 * the access point to the physics processes selection and customization. It
 * uses other classes (GeneralPhysics, EMPhysics, HadronPhysics, MuonPhysics,
 * IonPhysics) to setup particular physics processes. This part is quite recent
 * and has to be yet fully validated. Since the Geant4 toolkit is open to new
 * physics processes (along with new description of already present processes),
 * this will the access point for further development in the physics selection
 * (and thanks to good OO architecture of Geant4 this should be decoupled by
 * from the other part of the code)
 * - PrimaryGeneratorAction: this class, derived by a standard
 * G4VUserPrimaryGeneratorAction class, is in charge of production and iniection
 * of primary particles in the detector simulation. In our cases it is linked,
 * via its public interface, to the GAUDI service FulxSvc that is responsable to
 * generate the incoming fluxes of particles in an abstract and customizable
 * way. 
 * - DetectorManager: this class, that derives from a standard
 * G4VSensitiveDetector, manages the setup and working of the sensitive
 * detectors of the simulation and their interaction with the Transient Data
 * Store of GAUDI. In our case it is an abstract class, concretely implemented
 * in the two subclasses PosDetectorManager and IntDetectorManager; the first
 * one is associated with detector that saves hits information usign the
 * McPositionHit of GlastEvent (silicon planes in the TKR), while the second one
 * is used for detectors using the McIntegratingHit (ACD tiles and CAL cells).
 *
 * A special class DisplayManager uses the GuiSvc to display: 
 *          - All sensitive detectors, ACD tiles, TKR planes, CAL logs , 
 *            CAL diodes.  
 *          - Ids of sensitive detectors 
 *          - Hit detectors only 
 *          - Steps 
 *          - Tracks
 *
 * <br> The display of each can be enabled independently.  It will be replaced
 * or supplemented by HepRep/WIRED.  <br>
 *
 * Few other classes are used to help in the filling of the McParticle tree (the
 * TrackingAction and McParticleManager) and to deal with some oddity in the
 * output of Geant4 (UIsession)<br>
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
 * @todo Insert documentation for the new physics
 * @todo Insert new code for physics customization via jobOptions file
 *
 */

