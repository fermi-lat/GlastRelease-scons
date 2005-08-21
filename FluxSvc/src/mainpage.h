// $Header$
// Mainpage for doxygen

/*! \mainpage package FluxSvc

   \authors Toby Burnett, Sean Robinson, Theodore Hierath, and others.

 \section intro Introduction
  This package implements the Gaudi objects:
  
  - FluxSvc: mangages (using the package flux) the sources. Also implements the Runnable interface
  - FluxAlg: generates the current incoming particl
  - ExposureAlg: if a Clock (or "TimeTick") pseudo-particle was geneated, puts the info into the exposure history
  <br>
<br>
Usage is primarily via the FluxAlg algorithm, which access the service to generate a particle and place it in the TDS.
    
  <hr>
  @section fluxalg_jobOptions FluxAlg jobOptions
    @param FluxSvc.source_lib   [\$(FLUXSVCROOT)/xml/source_library.xml] list of file names containing source_lib elements
    @param FluxSvc.dtd_file     [\$(FLUXSVCROOT)/xml/source.dtd]  DTD file used to parse the XML files 
    @param FluxSvc.EvtMax       [0]  If non-zero, used as a maximum in the FluxSvc loop
    @param FluxSvc.StartTimeEnvVar [""] If set (e.g., to "runName") use this as the source of the start time.
    It has two forms: a number, or optionally a number for the start time offset, followed by a comma ant the delta time
    @param FluxSvc.StartTime    [0]  Mission elapsed time start
    @param FluxSvc.EndTime      [0]  Mission elapsed time end--if non-zero, will be the end.
    @param FluxSvc.DeltaTime    [0]  Maximum elapsed time.
    @param FluxSvc.StartDate    [""] Date, in form "2004-07-25 18:00" 
    note that if both StartDate and StartTime are specified, the latter is added to form the actual start

    @param FluxAlg.source_name  ["default"] source name, name must be in the source_lib files
    @param FluxAlg.sources      [{}] if used, can specify multiple sources. Overrides source_name
    @param FluxAlg.MCrun        [100] Initial run number
    @param FluxAlg.area         [6.0] target area in m^2
    @param FluxAlg.pointing_mode [0]  Corresponds to the following, from GPS
        - 0 No rocking rotation done at all.
        - 1 Satellite will be rocked toward the north pole in the northern hemisphere, opposite in the south.
        - 2 (experimental) like UPDOWN, except that rotation at equator happens gradually.
        - 3 LAT rocked northward for one orbit, southward for the next.
        - 4 fixed.
    @param FluxAlg.rocking_angle [0 deg] Rotation angle for Glast, about x-axis. 
    @param FluxAlg.rocking_angle_z [0 deg] Rotation angle for Glast, about z-axis.
    @param FluxAlg.pointing_info_tree_name ["MeritTuple"] If set, copy "Pt" values to it. See point_info for definitions.
    @param FluxAlg.save_pointing_info [false] Set true to save all entries in the pointing tuple. Normally saved by merit.

  @section exposurealg_jobOptions ExposureAlg jobOptions
    @param ExposureAlg.root_tree ["pointing_history"] name for the root tree to be filled if there are clock ticks
    @param ExposureAlg.pointing_history_input_file [""] name for file to give to astro::GPS do define an input pointing history.


    
    <hr>
  @section Basic_XML_Sources Sources
  This is a limited selection of sources. See the contents of the source library file for a complete list.
  @param default
  0.1 GeV gamma-rays coming from the vertical local direction.  Used for default tests.
  @param albedo_gamma
  Source that represents the Earth horizon albedo 
  @param albedo_electronpositron
  Source that represents the splash and re-entrant albedo electrons and positrons
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
  @param cosmic_muons
  special source that mimics non-interacting cosmics

 <hr>
 The current complete list is:
 @verbatim
 albedo_electronpositron
albedo_electronpositron_total
albedo_electronpositronavg_total
albedo_electronpositronavghi
albedo_electronpositronavglow
albedo_electronpositronhi
albedo_electronpositronlow
albedo_gamma
albedo_proton_avg
albedo_proton_max
all_gamma
backgndavgpdr
backgndmaxpdr
backgndmix
chime
chimeavg
chimemax
chimemin
cosmic_muons
crab-galactic
crab-pulsed-pointed
cremeavg
default
diffuse
diffuse-100mev
electron
electronavg
electronmax
galcenter
gamma_100_gev_normal
gamma_100_gev_uniform
gamma_100_mev_uniform
gamma_10_gev_uniform
gamma_1_gev_30deg
gamma_1_gev_60deg
gamma_1_gev_normal
gamma_1_gev_uniform
normal_gamma
surface_muons
timetick
timetick30s
vela
vertical_muons
vertical_surface_muons
@endverbatim

This is extracted from 
@verbinclude "source_library.xml"

  <br>
  <h2> Defining an external source </h2>
    See the interface definition IRegisterSource for information on how to link code external to this package.


  <hr>
  \section notes release notes
  release.notes
  \section requirements requirements
  \include  requirements
  <hr>
   <hr>
  \todo Complete and recalibrate the CompositeDiffuse structure

*/

