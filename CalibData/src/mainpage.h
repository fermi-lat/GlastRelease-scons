// Mainpage for doxygen

/** @mainpage package CalibData
  @author Joanne Bogart
  @section intro Introduction
  This package contains data model for the transient detector store:
  definitions of the classes and of hierarchy within the store and
  assignment of class ids.  

  The model is a 3-tiered hierarchy, but the only nodes with significant
  associated data are the leaf nodes at the bottom.  There are no
  internal references among the nodes.  The expectation is that clients
  will only want to access the leaf nodes, and that the client is 
  prepared to handle all the data associated with such a node.
 
  <ol>
   <li>Top node is the root; has no data attached </li>
   <li>Second-level nodes correspond to calibration types, such as 
   CAL_LightAsym, TKR_HotChan, etc. They have only a small data object
   (class CalibCLIDNode, defined in the CalibSvc package) used for
   bookkeeping.  Its only data is single field containing the class ID
   of its child nodes (all children have the samle class ID, corresponding
   to the calibration type).</li>
   <li>Third-level (leaf) nodes have a calibration data set attached.
   Different child nodes of the same second-level node correspond to
   different calibration flavors.   Most often a second-level node will 
   have a single child with flavor = "vanilla". </li>
  </ol>

  The CalibModelSvc class makes certain information internal to the
  CalibData package safely available to CalibDataSvc, which is in  the
  CalibSvc package. 

  @section requirements requirements
  @include requirements
  <hr>
  @section notes release.notes
  release.notes
  <hr> 
  @todo    Make individual tds classes for different calibration types
  @todo    Make CalibModel (CondModel ??) class to describe hierarchy
           in the TDS.
 */

