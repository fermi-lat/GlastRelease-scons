// $Header$
// Mainpage for doxygen

/*! \mainpage package FluxSvc

   \authors Toby Burnett, Sean Robinson, Theodore Hierath, and others.

 \section intro Introduction
  This package implements a Gaudi service, encapsulating the flux package. Its only
  function is to return an IFlux object, whose methods are implemented by the flux package.
  <br>
  the file /src/test/jobOptions.txt holds information used for the implementation of FluxSvc.
  it also holds the standard set of strings representing xml filenames, thus allowing multiple 
  xml files to be used. <a href="../FluxSvcDoc2.htm">More Documentation</a>
  <br>
  <h2> Defining an external source </h2>
    See the interface definition IRegisterSource.
    

  <hr>
  \section notes release notes
  release.notes
  \section requirements requirements
  \include requirements
  <hr> 
  \todo Complete and recalibrate the CompositeDiffuse structure
  \todo Overhaul time handling to use the TimeStamp class.

*/

