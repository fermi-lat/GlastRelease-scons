// $Header$
// Mainpage for doxygen

/** @mainpage package G4Propagator
 *
 * @authors Tracy Usher
 *
 * @section description Description
 *
 * This package is the GLAST interface to the Geant4 version of the "Particle
 * Propagator." The job of the propagator is to "propagate" the track parameters,
 * including error matrix, from the point at which they are defined to a given
 * end point. 
 *
 * Several classes are defined in order to accomplish this task:
 * 1) G4PropagatorTool - This is a Gaudi tool for interfacing external Glast packages
 *    to the propagator. It is also responsible for obtaining the necessary Geant4 
 *    geometry information from the G4GeoemetrySvc as well as obtaining a pointer to
 *    the Geant4 TrackingInformaton class
 *
 * 2) G4ParticlePropagator - The top level class which implements the IKalmanParticle
 *    interface for the propagator. Mostly implements the functions for returning specific
 *    values
 *
 * 3) ParticleTransporter - internal class which interfaces to the Geant4 volume code to 
 *    perform the actual transportation through the Geant4 geometry
 * 
 * 4) TransportStepInfo - an internal untility routine keeping track of information for
 *    each step.
 *
 *<hr>
 * @section jobOptions jobOptions
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
 * @todo Complete the implementation of "position"
 * @todo Determine the source of problem between TkrRecon and Propagator which causes crash
 * @todo Insert documentation
 *
 */

