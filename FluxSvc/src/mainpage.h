// $Header$
// Mainpage for doxygen

/*! \mainpage package FluxSvc

   \authors Toby Burnett, Sean Robinson, Theodore Hierath, and others.

 \section intro Introduction
  This package implements the Gaudi objects:
  
  - FluxSvc: mangages (using the package flux) the sources. Also implements the Runnable interface
  - FluxAlg: generates the current incoming particl
  - ExposureAlg: if a Clock (or "TimeTick") pseudo-particle was geneated, puts the info into the exposure history
  - PointInfoAlg; Manages calculation of position and orientation entries in the "merit" tuple, with an entry per
  event, and the pointing info tuple, fed by ExposureAlg.
  <br>
<br>
Usage is primarily via the FluxAlg algorithm, which access the service to generate a particle and place it in the TDS.
    
  <hr>
  @section fluxsvc_jobOptions FluxSvc jobOptions
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
    @param FluxSvc.SampleInterval [1.0] Governs the minimum time that must elapse before a the GPS observers
    are notified to update pointing and attitude information.
    @param FluxSvc.OrbitInclination [25.6] Orbit inclination angle, in degrees for default orbit
    @param FluxSvc.SAApolyLat  [{}] list of latitudes defining the SAA exclusion polygon
    @param FluxSvc.SAApolyLon  [{}] longitudes
    @param FluxSvc.xmlListFile    [""] name of a file containing a list of xml file names, multiple source_lib entries, consistent with obssim
    @param FluxSvc.EnableAberration [false] enable stellar aberration for celestial sources

  @section fluxalg_jobOptions FluxAlg jobOptions
    @param FluxAlg.source_name  ["default"] source name, name must be in the source_lib files
    @param FluxAlg.sources      [{}] if used, can specify multiple sources. Overrides source_name
    @param FluxAlg.MCrun        [100] Initial run number
    @param FluxAlg.area         [6.0] target area in m^2
    @param FluxAlg.backoff      [2.0] backoff distance in m
    @param FluxAlg.rocking_angle [0 deg] Rotation angle for Glast, about E-W axis. If set, so to 35 degress, and 
    there is no pointing history file, the mode will be to use the interal orbit description, with 
    rocking on alternate orbits. 
    @param FluxAlg.alignment    [{}] Set three rotation angles to align GLAST
    @param FluxAlg.misalignment [{}] Set three rotation angles to mis-align GLAST
    @param FluxAlg.pointingDirection [{}] Set (ra,dec) for pointed mode. Other rocking stuff ignored.
    @param FluxAlg.AvoidSAA     [false] set true to skip events during SAA interval
    @param FluxAlg.zenithTheta  [-99] if overridden, set to this angle in local zenith frame
    @param FluxAlg.PointingHistory [{}] Up to three strings: 
                                 (1) file name (text or FITS FT2),  
                                 (2) date-time offset needed for text format, 
                                 (3) any string, if present, will trigger the horizontal orientation
    @param FluxAlg.filterCone   [{}] Triplet: ra, dec, radius of cone to apply to celestial sources
    @param FluxAlg.sourceListFile [""] File with a list of source names to be added to the "sources" list. Consistent with obssim
    @param FluxAlg.abortOnExceptioni [false] What to do if flux generates an exception.

  @section exposurealg_jobOptions ExposureAlg jobOptions
    @param ExposureAlg.root_tree ["pointing_history"] name for the root tree to be filled if there are clock ticks
    @param ExposureAlg.PrintFrequency [1] number of INFO entries to create per tick
    @param ExposureAlg.clockName ["clock"]  The source name that will trigger entries

  @section PointInfoAlg_jobOptions PointInfoAlg jobOptions
    @param PointInfoAlg.pointing_info_tree_name ["MeritTuple"]If set, copy "Pt" values to it. See point_info for definitions.
    aparam PointInfoAlg.save_pointing_info  [false] Set true to save all entries in the pointing tuple. Normally saved by merit.

  @section OrbitSvc_jobOptions OrbitSvc jobOptions
    @param OrbitSvc.zenithTheta  [-99] if overridden, set to this angle in local zenith frame
    @param OrbitSvc.pointingDirection [{}] Set (ra,dec) for pointed mode. Other rocking stuff ignored.
    @param OrbitSvc.PointingHistory [{}] Up to three strings: 
                                 (1) file name (text or FITS FT2),  
                                 (2) date-time offset needed for text format, 
                                 (3) any string, if present, will trigger the horizontal orientation
    @param OrbitSvc.rocking_angle [0 deg] Rotation angle for Glast, about E-W axis. If set, say to 35 degress, and 
    there is no pointing history file, the mode will be to use the interal orbit description, with 
    rocking on alternate orbits. 

    <hr>
  @section Basic_XML_Sources Sources

  The sources names and definitions are defined by the XML files that are loaded. See xml/source_library.xml in the flux
  package, and xml/source_library.xml in the CRflux package. 
 <hr>

@endverbatim

This is extracted from 
@verbinclude "xml/source_library.xml"

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

