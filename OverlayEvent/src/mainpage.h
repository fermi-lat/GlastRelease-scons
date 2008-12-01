/** @mainpage  
 package OverlayEvent

 @section Introduction Introduction
 This package contains definitions for all GLAST Overlay data to be stored on 
 the Transient Data Store (TDS). Currently, the OverlayEventModel defines that
 these objects be stored in the "/Event" section of the TDS which has the 
 advantage of controlling the clearing/destruction of the objects automatically
 at begin event time. 

 ******** NOTE: this will change soon *************

 Note that all inherit from DataObject, and correspond to transient 
 store objects at the top level, under "/Event".

 @section EventOverlay EventOverlay

  This is a top level data object for Overlay that contains the identical information
  as the Event package's Event object, but the EventOverlay object will contain the 
  run/event/etc. information specific to the input background event, not the current
  simulated event.

 @section AcdOverlay AcdOverlay
  
  Overlay data for the ACD

 @section CalOverlay CalOverlay
  
  Overlay data for the CAL

 @section TkrOverlay TkrOverlay
  
  Overlay data for the TKR

 @section GemOverlay GemOverlay
  
  Overlay data for the GEM

 <hr>
 @section jobOptions jobOptions
 No jobOptions are used within this package.
 <hr>
 @section Tests Tests
 There are no test routines in this package:
 <hr>
 @section notes release notes
 release.notes
 <hr>
 @section requirements requirements
  @verbinclude requirements
*/
