/** @mainpage  
* package Event
*
* @section Introduction Introduction
* This package contains definitions for all GLAST event data to be stored on 
* the TDS:
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
*   -# CalXtalRec Collection
*   -# CalCluster Collection
* - TkrRecon
*   -# TkrCluster Collection
*   -# TkrPatCand Collection
*   -# TkrFitTrack Collection
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
* The Monte Carlo data is stored on the following path on the TDS:
* /Event/MC/
*
* There are four Monte Carlo TDS classes:
* - McParticle
* - McPositionHit
* - McIntegratingHit
* - McTkrStrip
*
* @section DigiEvent DigiEvent
* The digitization classes represent our detector data.  This data appears on the
* TDS under /Event/Digi/
*
* There are three Digi TDS classes:
* - AcdDigi
* - CalDigi
* - TkrDigi
*
* @section Reconstruction Reconstruction
* The reconstruction data is available on the TDS.  This data appears under a 
* variety of headings depending upon the subsystem:
* /Event/AcdRecon
* /Event/CalRecon
* /Event/TkrRecon
*
* @section reference Reference Document
*  See the formal  
*  <a href="http://www.slac.stanford.edu/~hansl/glast/note/rd.pdf">documentation</a>.
*
*  <hr>
* @section notes release notes
* release.notes
* <hr>
* @section requirements requirements
* @include requirements
*/
