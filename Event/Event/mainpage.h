//$Header$
// (Special "header" just for doxygen)

/*! @mainpage  package GlastEvent

This package contains definitions for all GLAST event data:
	- Event
	- MCEvent
	- RawEvent
	- RecEvent
	- Run

Note that all inherit from DataObject, and correspond to transient 
store objects at the top level, under "Event".

\section MCEvent MCEvent
The MCevent corresponds to all the Hit classes, and the MC truth. 
- DataObject classes
  - CalorimeterHits
  - GlastHits
  - SiLayerHits
  - TowerHits
  - TrackerHits
- ContainedObject classes
  - ACDTileHits
  - SiStripHits
  - CalorimeterLogHits

All of the above are "wrapper" classes, inheriting from GlastEvent::HitsBase. They are 
instantiated with concrete classes from the <b>instrument</b> package. See the Gsim*.h headers:
  -GsimACDTileHits.h
  -GsimCalorimeterHits.h
  -GsimCalorimeterLogHits.h
  -GsimGlastHits.h
  -GsimSiLayerHits.h
  -GsimSiStripHits.h
  -GsimTowerHits.h
  -GsimTrackerHits.h

  An exception is the GlastEvent::EventHits class, inheriting from DataObject, and containing
  GlastEvent::GlastHits.

\section reference Reference Document
  See the formal  
  <a href="http://www.slac.stanford.edu/~hansl/glast/note/rd.pdf">documentation</a>.

  <hr>
  \section notes release notes
  \include release.notes
  \section requirements requirements
  \include requirements


*/
