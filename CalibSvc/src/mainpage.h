// Mainpage for doxygen

/** @mainpage package CalibSvc
  @author Joanne Bogart
  @section intro Introduction
  This package contains classes and interfaces for calibration-related
  services, in particular a data service, a conversion service for
  the metadata, and (ultimately) conversion services with associated
  converters for individual calibration classes.

  CalibDataSvc inherits from the Gaudi class DataSvc.  Value added
  includes implementation of Gaudi interface IDetDataSvc 
  (maintenance of event time, so that particular calibration data sets 
  may be compared with it to check validity for the current event) and of
  abstract interface IInstrumentName, defined in this package, which
  plays a similar role but for instrument rather than time.

  CalibDataSvc::initialize() sets up the non-leaf nodes for the calibration
  TDS and registers addresses (but not objects) for the leaf nodes.

  CalibMySQLCnvSvc is a conversion service for the metadata (i.e., satisfies
  the ICalibMetaCnvSvc abstract interface defined in this package) when
  the metadata has persistent form as a row in a MySQL table.  The service
  has no associated converters.  It overrides the standard conversion service
  handling of converter methods (which is to delegate to a converter) and
  handles these methods itself.  It does not register any data objects
  in the calibration TDS.  Instead, its createObj method

    <ul>
      <li>Finds the best-match row in the MySQL database for the calibration
          requested</li>
      <li>Fetches information from that row needed to access the bulk data
          for the calibration, namely fields describing persistent format
          (ROOT or XML), file identifier, and format version.  It also
          fetches other information such as start and stop of validity
          interval for this calibration.</li>
      <li>Forms a new opaque address, using this information, and 
          invokes the persistency service createObj method with this
          new address.  This should cause the proper conversion service
          (ROOT or XML) to be invoked, which should in turn delegate
          conversion to an appropriate converter, which will form and
          register the object.</li>
    </ul> 

  
  The package will also contain CalibROOTCnvSvc and CalibXMLCnvSvc
  classes and their converters, to be written.

  @section requirements requirements
  @include requirements
  <hr>
  @section notes release.notes
  release.notes
  <hr> 

  @todo    Figure out what other information from the metadata needs to
           be acquired and saved, and where it should go.
  @todo    Write CalibROOTCnvSvc and CalibXMLCnvSvc
  @todo    Define calibration TDS classes in detail; write converters.
           
 */

