/** @mainpage  
* package Event
*
* @section Introduction Introduction
* This package contains definitions for all GLAST event data to be stored on 
* the Transient Data Store.  During each event, data is created and shared on
* the TDS.  At the end of each event, this data is cleared from the data store.
* The TDS data currently includes:
* - Event
* - MCEvent
*   -# McParticle Collection
*   -# McPositionHit Collection
*   -# McIntegratingHit Collection
* - DigiEvent
*   -# AcdDigi Collection
*   -# CalDigi Collection
*   -# TkrDigi Collection
* - AcdRecon
* - CalRecon
*   -# CalXtalRecData Collection
*   -# CalCluster Collection
* - TkrRecon
*   -# TkrCluster Collection
*   -# TkrPatCand Collection
*   -# TkrTrack Collection
*   -# TkrVertex Collection
*
* Note that all inherit from DataObject, and correspond to transient 
* store objects at the top level, under "/Event".
*
* A description of the GLAST Event Model is available on the web:
* <A HREF="http://www-glast.slac.stanford.edu/software/gaudi/tds/event_model.htm> GLAST Event Model </A>
*
* @section Event Event
* The top-level data object on the TDS contains the Event header information,
* stored in a class called EventHeader.  The data includes:
* - Run Number
* - Event Id
* - Time Stamp
* - Trigger Word
*
* @section MonteCarlo MonteCarlo
* The Monte Carlo data is stored on the TDS path:
* /Event/MC/
*
* The top-level Monte Carlo header class is called MCEvent.  This class contains
* data member, a source id.  
*
* Under the MC header, there are four Monte Carlo TDS classes:
* - McParticle Collection
* - McPositionHit Collection
* - McIntegratingHit Collection
* - McTkrStrip Collection
*
* @section DigiEvent DigiEvent
* The digitization classes represent our detector data.  This data appears on the
* TDS under /Event/Digi/
* There is one data member in the Digi header, a flag denoting whether or not 
* digitization data origined from Monte Carlo data or not.
*
* Under the Digi header, there are three Digi TDS classes:
* - AcdDigi Collection
* - CalDigi Collection
* - TkrDigi Collection
*
* @section Reconstruction Reconstruction
* The reconstruction data is available on the TDS.  This data appears under a 
* variety of headings depending upon the subsystem:
* - /Event/AcdRecon
* - /Event/CalRecon
*    -# CalXtalRecData Collection
*    -# CalCluster Collection
* - /Event/TkrRecon
*    -# TkrCluster Collection
*    -# TkrPatCand Collection
*    -# TkrTrack Collection
*    -# TkrVertex Collection
*
* @section reference Reference Document
*  See the formal  
*  <a href="http://www.slac.stanford.edu/~hansl/glast/note/rd.pdf">documentation</a>.
*
* @section jobOptions jobOptions
* No jobOptions are used within this package.
*
* @section Tests Tests
* There are two test routines in this package:
* - test_Event
* - test_Tables
*
*  <hr>
* @section notes release notes
* release.notes
* <hr>
* @section requirements requirements
* @verbinclude requirements
*/
