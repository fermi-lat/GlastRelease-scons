// Mainpage for doxygen

/** @mainpage package calibUtil
 
  @authors J. Bogart, M. Gulati
 
  @section intro Introduction
  This package contains (well, will contain) utilities for storing,
  searching for, and accessing calibration data.  "Calibration data"
  include at a minimum hardware status information, such as hot and
  dead tracker strips, and conversion constants, such as those 
  characterizing light attenuation in calorimeter crystals.
 
  The overall strategy is to store a limited amount of <b>metadata</b> 
  about each set of calibration data in a relational database table, 
  while the data are kept in flat files, one per result set.  Most 
  likely two formats will be supported for these files: XML and ROOT.
  See the <a href="http://www.slac.stanford.edu/exp/glast/ground/software/calibration/CalibSvcSpec.shtml">Calibration Services Specification</a> document
  for more details.
 
  The classes in this package, all in the calibUtil namespace, include
  
  <dl>
   <dt>Metadata</dt>
   <dd>Provides services to read and write metadata.  All access to the
       rdbms (currently MySQL) is encapsulated here.</dd>
   <dt>StripSrv</dt>
   <dd>Reads in bad (tracker) strip data from its persistent form (currently
       XML) and makes it available to clients via query services and
       a visitor serice, traverseInfo. </dd>
   <dt>ClientObject</dt>
   <dd>This pure virtual class is used by the StripSrv visitor.  Clients
       must implement their own derived class.</dd>
   <dt>GenericSrv</dt>
   <dd>Handles generic part of calibration data set; that is, that part
       which is common to all types of calibration.  This is an internal
       calibUtil class not directly accessible by calibUtil clients.</dd>
   <dt>Timestamp</dt>
   <dd>Encapsulation of timestamp (resolution of a second) as used within
       calibUtil.</dd>
   </dl>
  
 
  <hr>
  @section notes release.notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements
  <hr> 
  @todo 
     Finish Metadata searching services

  @todo 
     Add write services to StripSrv

  @todo 
     Make analagous class(es) for handling ROOT persistent form
     
 */

