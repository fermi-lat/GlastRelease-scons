
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
* The CLHEP Random engine to use - same for all Dll's that use random
* numbers.
* Default:TripleRand
* @param GlastRandomSvc.SeedFile
* An ASCII file containing a set of lines containing
* space separated run and sequence numbers to use. If this
* file is not in the same directory as Gleam.exe it's full
* path must be specified. Note that specification of this
* file overrides individual run and initial sequence numbers
* set via GlastRandomSvc.RunNumber and 
* GlastRandomSvc.InitialSequenceNumber.  
* @param GlastRandomSvc.RunNumber
* Run number. Note that if a file of run and sequence numbers is 
* specified via the job options parameter GlastRandomSvc.SeedFile 
* that overrides individually set run and initial sequence 
* numbers.
* Default:10
* @param GlastRandomSvc.InitialSequenceNumber
* Incident particle sequence number to begin run with - is not
* event number as not all events trigger.
* Note that if a file of run and sequence numbers is specified 
* via the job options parameter GlastRandomSvc.SeedFile 
* that overrides individually set run and initial sequence 
* numbers.
* Default:0
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

*/

