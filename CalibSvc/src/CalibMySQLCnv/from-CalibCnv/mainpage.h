// Mainpage for doxygen

/** @mainpage package CalibData
  @author Joanne Bogart
  @section intro Introduction
  This package contains data model for the transient detector store:
  definitions of the classes and of hierarchy within the store.

  @section requirements requirements
  @include requirements
  <hr>
  @section notes release.notes
  release.notes
  <hr> 
  @todo   Write persistency service class (DetCnvSvc ??)
  @todo   Write actual converters.  Probably there is just one top-level
          one for calibration data, which finds the "right" record in 
          the metadata database; then invokes something else to convert 
          the persistent form of the data and register it in the tds
 */

