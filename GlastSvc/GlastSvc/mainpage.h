// $Header$
// Mainpage for doxygen

/*! \mainpage package GlastSvc

\section NewConverters New Converters

Heather and Ian have defined a new set of converters to go from the IRF format to
the new set of Monte Carlo classes. They all inherit from from the BaseCnv class.

  - BaseCnv classes
    - MCACDHitCnv
    - MCCalorimeterHitCnv
    - MCSiLayerCnv

\section Testing Testing

Heather has made a small testing Algorithm called CreateEvent.

\section GlastEventConverter Changes to GlastEvent Converter

  The GlastEventSelector was modified (02/28/2001) to handle the case of Input set
  to NONE. i.e.

  EventSelector.Input = "NONE" is now ok.

  Also a basic counter was added to GlastEvtIterator that enable the event loop to 
  be terminated from the GUI.

Also working in conjuction with the MC Converters is the IRFConverter which uses 
the IRF parsers defined in the instrument package.

\section General General


 Defines the services
 - IGlastDetSvc, implmented by GlastDetSvc. 
 - IGlastIRFLoadSvc also implemented by GlastDetSvc.

  <hr>
  \section notes release notes
  \include release.notes
  \section requirements requirements
  \include requirements

*/

