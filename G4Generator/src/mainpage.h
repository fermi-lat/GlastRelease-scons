// $Header$
// Mainpage for doxygen

/** @mainpage package G4Generator
 *
 * @authors T.Burnett, R.Giannitrapani, F. Longo
 *
 * @section description Description
 *
 * This package is the GLAST interface to the Monte Carlo simulations toolkit
 * Geant4 (see <a href="http://geant4.web.cern.ch/geant4/">here</a> 
 * for more information on this C++ toolkit). 
 * The steering is done with the Gaudi algorithm G4Generator.<br>
 * 
 * The main simulation is controlled by RunManager, a customized and stripped
 * version of the Geant4 standard G4RunManager singleton. Since GLAST main event
 * loop is driven by GAUDI and it will not use any graphics or data persistency
 * features of Geant4, we have included in the RunManager only the real
 * necessary parts for setup and run the generator. RunManager itself uses the
 * following classes
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
 * this will be the access point for further development in the physics
 * selection (and thanks to good OO architecture of Geant4 this should be
 * decoupled from the other part of the code)
 * - PrimaryGeneratorAction: this class, derived by a standard
 * G4VUserPrimaryGeneratorAction class, is in charge of production and iniection
 * of primary particles in the detector simulation. In our cases it is linked,
 * via its public interface, to the GAUDI service FulxSvc that is responsable to
 * generate the incoming fluxes of particles in an abstract and customizable
 * way. A new design is on the way by which the primary particle will be created
 * by an algorithm (FluxAlg) using FluxSvc and stored in the TDS as the root of
 * an McParticle tree; G4Generator will than retrive this primary particle from
 * the TDS and initialize the PrimaryGeneratorAction with that.
 * - DetectorManager: this class, that derives from a standard
 * G4VSensitiveDetector, manages the setup and working of the sensitive
 * detectors of the simulation and their interaction with the Transient Data
 * Store of GAUDI. In our case it is an abstract class, concretely implemented
 * in the two subclasses PosDetectorManager and IntDetectorManager; the first
 * one is associated with detector that saves hits information usign the
 * McPositionHit of GlastEvent (silicon planes in the TKR), while the second one
 * is used for detectors using the McIntegratingHit (ACD tiles and CAL cells).
 *
 * A special class DisplayManager uses the GuiSvc to display (optionally): 
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
 * TrackingAction and McParticleManager), to manage the Geant4 transport of
 * particles for reconstruction uses (G4PropagatorSvc, TransportStepInfo,
 * ParticleTransporter and G4ParticlePropagator) and to deal with some oddity in
 * the output of Geant4 (UIsession)<br>
 * 
 * @section pruning McParticle tree pruning
 *
 * Since some events has an McParticle tree too big to be easily browsed and
 * anlyzed, it has been necessary to design some sort of pruning mechanism for
 * such a tree; this pruning is performed at the end of the event and affect the
 * TDS content. When an McParticle is pruned, a special care must be taken for
 * the parentage relationships; in our case the McParticle pruned is removed
 * from the daughters of its mother, while all the daughters of the pruned
 * McParticle are added to the daughters list of its mother.
 *
 * The parameter <i>G4Generator.mcTreeMode</i> can be used to set the level of
 * pruning of the McParticle tree; there are three possible values:
 *
 * - <b>full</b> in this case the full McParticle is saved in the TDS.
 * - <b>pruneCal</b> in this case all the particles that are created and end in a
 *                   volume belonging to the CAL are pruned, unless they generate a PosHit. 
 * - <b>minimal</b> in this case only the primary is retained in the tree and
 *                  all the particles "generated" by it.
 *
 *<hr>
 * @section jobOptions jobOptions
 * @param G4Generator.UICommands  default value ""
 * @param G4Generator.geometryMode [""] Blank means use default from GlastDetSvc
 * @param G4Generator.saveTrajectories [0] Set non-zero to allow trajectory objects (Event::McTrajectory) to be saved on the TDS
 * @param G4Generator.mcTreeMode ["minimal"]
 * @param G4Generator.defaultCutValue [0.1*mm] The cut to be applied to all physics processes, all materials
 * @param G4Generator.physics_choice ["full"]
 * @param G4Generator.mscatOption [true]  If true: replace the native G4 multiple scattering code with a 
 different version managed by Geant4/MultipleScatteringFactory.
 This can use 
any particle (apart only I think the lambda test ..) that one has from 
FluxSvc. Specyfying only muons, I think you could get a segmentation fault 
in the case you will use the proton for example. The "only-em" has the Decay processes, the muon 
processes (where I commented out the MuonAbsorption), all the e-/e+ and 
gamma processes and the multiple scattering + ionisation for hadrons and 
ions. For the full there is the actual choice of processes. 
 * @param G4Generator.eLossCurrent [true]  If true: replace the native G4 energy loss code with the 5.2 version.


 *
 * @section test test
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

