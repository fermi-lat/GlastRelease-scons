/** @mainpage package GlastSvc
*
* @section Introduction
* This package defines GLAST-specific Gaudi services and interfaces.
*
* @section EventCnvSvc EventCnvSvc
* This is GLAST's Persistency Service.  
* A Persistency Service requires:
* - A service derived from IConversionSvc
* - A set of converters
* - IOpaqueAddresses to associate converters with paths on the TDS
*
* For information about converters and the persistency service see the Gaudi
* Developers Guide <A HREF="http://proj-gaudi.web.cern.ch/proj-gaudi/GDG/v2/Output/GDG_Converters.html#1010951"> Chapter 13</A>
*
* At this point, the service does not do very much except setup simple 
* converters for the Event and McEvent header objects.  Currently, we are not
* using the persistency service to handle our persistent store.  This will 
* be added in the future.
*
* The set of converters available to the EventCnvSvc are all derived from
* the BaseCnv class.  Currently, there are 2 concrete converters:
* - EventCnv
* - MCEventCnv
*
* @section GlastDetSvc GlastDetSvc
* A GLAST-defined Gaudi service which provides access to detModel through a 
* wrapper class, DMmanager.  GlastDetSvc thus provides access to geometry
* parameters and constants.
*
* @section GlastRandomSvc GlastRandomSvc
* A GLAST-defined Gaudi service which provides for the setting of CLHEP
* engines on a per Dll basis and provides for the setting of seeds
* for the random number generators for each Dll on a per incident
* particle basis.
*
* @section CurrentEventAlg CurrentEventAlg
* An algorithm which creates an output ASCII file and writes the run and event
* ids for the event currently being processed.  The file is flushed for each
* event, in case a crash occurs during processing.
*
* @section Interfaces Interfaces
* A number of interfaces are defined in this package:
* - IEventCnvSvc
* - IGeometry
* - IGlastDetSvc
* - IKalmanParticle
* - IMedia
* - IPropagatorSvc
* - IRecoSvc
*
* @section jobOptions jobOptions
* @param GlastDetSvc.xmlfile
* Name of the top-level XML file that contains the geometry description
* @param GlastDetSvc.topVolume
* Top-level volume to create for this run  This allow one to define
* a subset of the GLAST instrument for study.
* Default:  LAT
* @param GlastDetSvc.visitorMode
* Default: propagate
* @param GlastRandomSvc.RandomEngine
* The CLHEP Random engine to use - same for all Dll's that use random numbers.
* Default:TripleRand
* @param GlastRandomSvc.SeedFile
* Currently, this parameter has no effect!
* The name SeedFile is misleading anyway!
* An ASCII file containing a set of lines containing space separated run and
* sequence numbers to use. If this file is not in the same directory as
* Gleam.exe it's full path must be specified. Note that specification of this
* file overrides individual run and initial sequence numbers set via
* GlastRandomSvc.RunNumber and GlastRandomSvc.InitialSequenceNumber.
* @param GlastRandomSvc.Seed
* This parameter should be used in emergency situations only!
* The seed to be used. In the case of negative and/or even values a warning is
* logged. Overwrites the seed calculated by GlastRandomSvc::applySeeds() if all
* following conditions are true:
* *) GlastRandomSvc::autoSeed=false, i.e. run seeded runs only
* *) one random engine instantiation only, i.e. valid for Linux/gcc only
* *) GlastRandomSvc.Seed!=0
* Default:0
* @param GlastRandomSvc.RunNumber
* Run number.
* Note that if a file of run and sequence numbers is specified via
* GlastRandomSvc.SeedFile that overrides individually set run and initial
* sequence numbers.
* Default:10
* @param GlastRandomSvc.RunNumberString
* Run number in string format, allowing the use of environment variables. If set
* non-null, it takes precedence over the integer version.
* Note that if a file of run and sequence numbers is specified via
* GlastRandomSvc.SeedFile that overrides individually set run and initial
* sequence numbers.
* Default:null
* @param GlastRandomSvc.InitialSequenceNumber
* Incident particle sequence number to begin run with - is not event number as
* not all events trigger.
* Note that if a file of run and sequence numbers is specified via
* GlastRandomSvc.SeedFile that overrides individually set run and initial
* sequence numbers.
* Default:0
* @param GlastRandomSvc.NumberOfRuns
* Specifies the number of runs contained in a task. This number is being used in
* the algorithm which calculates the seed, and ensures that all events of all
* runs are distinct.  All runs of a task have to run with the same
* GlastRandomSvc.NumberOfRuns!
* Default:20000 for  GlastRandomSvc.autoSeed=true, 20000000 (20M) for false
* @param GlastRandomSvc.autoSeed
* A better name would be eventSeeded.
* If set to true, GlastRandomSvc::applySeeds() is called for each event, using
* the run number,  GlastRandomSvc.NumberOfRuns, and the current event id. Thus,
* each event generation starts with it's own well defined seed. Using odd
* positive seeds only, this allows for 2^30 (1Gi) distinct events.
* If set to false, GlastRandomSvc::applySeeds() is called only once at the
* beginning of a run, using the run number, GlastRandomSvc.NumberOfRuns, and
* GlastRandomSvc.InitialSequenceNumber. Using odd positive seeds only, this
* allows for 1Gi distinct runs, each with up to 10^42 events (the latter was
* calculated assuming 100k random numbers for each event).
* Default:false
* @param GlastRandomSvc.pedantic
* GlastRandomSvc checks various parameters for validity:
* *) the run number has to be positive smaller than GlastRandomSvc.NumberOfRuns
* *) GlastRandomSvc.InitialSequenceNumber has to be positive
* *) in the case of GlastRandomSvc.autoSeed=true, checks if for given
*    GlastRandomSvc.NumberOfRuns the id of the last event stays within limits to
*    avoid duplicate events. I.e., the following inequality has to be valid:
*    GlastRandomSvc.NumberOfRuns*(lastEventId+1)<2^30.
* A failure in any of these checks causes a warning.  Additionally, if pedantic
* is set to true, GlastRandomSvc terminates with StatusCode::FAILURE.
* Default:true
* @param CurrentEventAlg.fileName
* Name of the output ASCII file which will contain the run and event ids 
* written by the CurrentEventAlg algorithm.  This JO parameter can handle
* environment variables as part of the input.  By default this fileName is "" 
* causing no eventId.txt file to be written.
*
* @section Testing test_GlastSvc
* The test routine contains one main algorithm called CreateEvent.
* Currently it tests access to GlastDetSvc during CreateEvent::initialize( ).
*
* <hr>
* @section notes release notes
* release.notes
* <hr>
* @section requirements requirements
* @verbinclude requirements
*
*/
