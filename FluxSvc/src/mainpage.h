// $Header$
// Mainpage for doxygen

/*! \mainpage package FluxSvc

   \authors Toby Burnett, Sean Robinson, Theodore Hierath, and others.

 \section intro Introduction
  This package implements a Gaudi service, encapsulating the flux package. Its only
  function is to return an IFlux object, whose methods are implemented by the flux package.
  <br>

  The FluxSvc package is designed to hold information on : 
  the orbital characteristics of GLAST,
  the Emissive characteristics of various sources, and
  the 'Current' time and a model for particle emission from desired sources.

  Through implementing a Gaudi Service (FluxSvc.h), it:
  Takes user specifications on desired source/sources,
  Provides a user interface to produce additional incoming particles, and
  Provides the current time, particle energy, direction, and type of particle. 

  <br>
  the file /src/test/jobOptions.txt holds information used for the implementation of FluxSvc.
  it also holds the standard set of strings representing xml filenames, thus allowing multiple 
  xml files to be used.
  <a href="http://d0.phys.washington.edu/~srobinsn/FluxSvc/">See here for an explanation of how Gaudi was introduced to flux to make FluxSvc, and further documentation:</a>
  <br>
  <h2> Defining an external source </h2>
    See the interface definition IRegisterSource.
    

  <hr>
  \section notes release notes
  release.notes
  \section requirements requirements
  \include requirements
  <hr>
  @section jobOptions jobOptions
  @param FluxAlg.source_name
  This is the source name (defined in the XML) that is to be used as
  a source in the simulation.
  @param FluxSvc.source_lib
  this is a list of the XML files that will make their declared sources available to the 
  simulation.
  <hr>
  @section Basic_XML_Sources Basic_XML_Sources
  @param default
  0.1 GeV gamma-rays coming from the vertical local direction.  Used for default tests.
  @param albedo_gamma
  Source that represents the Earth horizon albedo with Glast zenith pointing
  @param albedo_electronpositron
  Source that represents the spalsh and re-entrant albedo electrons and positrons
  @param diffuse
  diffuse extraglactic from 10 MeV: from APJ 494:523
  @param diffuse-100mev
  Diffuse extraglactic from 100 MeV
  @param crab-galactic
  the Crab, pulsed portion, with pointed observation, for photons above 100 MeV,
  based on Nolan APJ 409:697
  @param electron
  galactic electron spectrum
  @param normal_gamma
  E^-1 spectrum from 18 MeV to 18 GeV and normal incidence
  @param muons
  special source that mimics non-interacting cosmics

  <hr>
  \todo Complete and recalibrate the CompositeDiffuse structure

*/

