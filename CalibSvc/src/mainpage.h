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

    <ol>
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
    </ol> 

  Several different modes of operation regarding use of event timestamps
  to select calibrations are supported:

    <ul>
      <li>Use event time from event; calibration validity interval must
          include this time.  Currently implemented only for data read out
          from a physical instrument, but support Monte Carlo data 
          expected soon.  </li>
      <li>Generate fake event time for each event and compare to validity
          interval as above. </li>
      <li>Don't use event time or validity interval at all.  Base calibration
          selection on the time at which the calibration was entered into
          the database.</li>
     </ul>

  See job options description below to learn how to select one of these
  modes.

  
  The package also contains CalibXMLCnvSvc and several converters for
  calibration data of different types stored in XML files.  These include
  (Tracker) hot and dead strips, Calorimeter gains, pedestals, integral
  non-linearity, and so forth.  CalibROOTCnvSvc and its converters are TBW.

  @section requirements requirements
  @include requirements
  <hr>
  @section notes release.notes
  release.notes
  <hr> 
  @section jobOptions jobOptions

  CalibDataSvc has the following job options properties:
  <dl>
  <dt> CalibStorageType</dt> <dd> defaults to MYSQL_StorageType</dd>
  <dt> CalibFlavorList</dt> <dd>list of additional flavors (beyond vanilla,
       which is always used) for which nodes should be created for each
       calibration type.   </dd> 
  <dt> CalibNameList</dt> <dd>list of calibration type/flavor
       combinations beyond those of flavor "vanilla" or of a flavor 
       specified in CalibFlavorList for which the data service is requested
       to create a node.</dd>
  <dt> CalibRootName</dt>  <dd>defaults to "Calib", top node in TCDS</dd>
  <dt> CalibInstrumentName</dt>  <dd>defaults to "LAT"</dd>
  <dt> UseEventTime</dt>         <dd>defaults to "true", corresponding
       to modes 1 and 2 above. If set to "false", 
       must also set CalibMySQLCnvSvc.UseEventTime to "false"</dd>
  <dt> CalibTimeSource</dt> <dd> Use value "data" for actual instrument
       data, "mc" (when implemented) for Monte Carlo, "clock" for fake
       event time.  Default is "none", also the correct value for mode 3.
       Currently it is not strictly necessary to specify "clock" for
       fake event time, but this implementation is expected to change.</dd>
  <dt> DbName</dt>         <dd>defaults to "calib", the production dbs for
       calibration metadata.  Algorithm developers, etc., may need to
       use the development database, "calib_user", instead.
  </dd>
  </dl>

  The service CalibMySQLCnvSvc has the following job options properties:
  <dl>
  <dt> Host </dt>        <dd>defaults to "*", meaning "use default MySQL
                             host" </dd>
  <dt> UseEventTime</dt>         <dd>defaults to "true". If set to "false", 
       must also set CalibDataSvc.UseEventTime to "false". In this case,
       calibrations will be selected according to enter_time (see further
       job options) rather than by validity interval compared with event time.
       </dd>
  <dt> EnterTimeStart</dt> <dd> Lower bound on calibration enter_time
Ignored unless UseEventTime is false.  </dd>
  <dt> EnterTimeEnd</dt> <dd> Upper bound on calibration enter_time.
       Default is "no upper bound". Ignored unless UseEventTime 
       is false.  </dd>
  </dl>


  The algorithm CalibEvtClock, used in the test program to generate fake event
  times and store them with CalibDataSvc, may also be scheduled by user
  applications.  It has job options properties

  <ul>
  <li>  startTime, defaults to current time. Time assigned to first event </li>
  <li>  delayTime, defaults to 2 seconds. Difference in timestamps between
        adjacent events</li>
  </ul>

  @todo    Figure out what other information from the metadata needs to
           be acquired and saved, and where it should go.
  @todo    Write CalibROOTCnvSvc
  @todo    Define remaining calibration TDS classes in detail; 
           write converters.
  @todo    implement CalibDataSvc.CalibTimeSource = "mc" 
  @todo    do fake event clock implementation as "data" is done (and "mc"
           will be done), without requiring user to configure a 
           separate algorithm in job options
 */

