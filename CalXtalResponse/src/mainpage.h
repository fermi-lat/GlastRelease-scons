// Mainpage for doxygen

/** @mainpage package CalXtalResponse
 *
 * @section CalCalibSvc CalCalibSvc
 * 
 *
 * ICalCalibSvc provides a Cal specific streamlined interface to the GLAST calorimeter calibration
 * database.  CalCalibSvc is currently the only concrete implementation of this interface.  
 * CalCalibSvc supports the following features:
 * 
 * - Simple interface: requires only specification of unique Cal xtal and range as input.  
 * Values are returned for the most part as C primitives (float, int, etc..) C++ std::vectors are used where appropriate.
 * -  Gleam/Gaudi mechanics such as TDS access, validity checking, calibration data format, database access, and data storage are all transparent to the user.
 * - Some calibration types are vectors which represent the 'knots' on a spline curve.  
 * CalCalibSvc generates the spline objects for these types for the user.  CalCalibSvc handles storage/deallocation/validation for these objects.  
 * Once created, spline objects are retained in memory until they become invalid.  At this point, the cache is flushed and new objects are created.
 * - CalCalibSvc can be configured to assign different flavors to each calibration type.
 * - Multiple instances of CalCalibSvc may be created if the user needs more than one flavor for one or more calibration types.
 *
 * @subsection jobOptions jobOptions
 * 
 * Suitable default values are currently provided for all options, all parameters are 'optional'.
 *
 * @param CalibDataSvc
 * Name of service from which calibration data is obtained. (default="CalibDataSvc")
 * @param DefaultFlavor
 * Flavor for all calibration types unless otherwise specified (default="vanilla")
 * @param FlavorIntNonlin
 * calib flavor for integral nonlinearity constants.
 * @param FlavorAsym
 * override calib flavor for asymmetry data.
 * @param FlavorPed
 * override calib flavor for pedestal data.
 * @param FlavorMeVPerDac
 * override calib flavor for MeVPerDac data.
 * @param FlavorTholdCI
 * override calib flavor for tholdCI data.
 * @param FlavorTholdMuon
 * override calib flavor for tholdMuon data.
 *
 * @section CalResponseTools CalResponseTools
 * 
 * The CalReponseTools are a collection of tools which provide a consistent interface to calorimeter response calculations.  
 * Once again, as much Gleam/Gaudi mechanics as possible are hidden from the user.
 *
 * CalResponseTools consists of 3 interfaces.
 * -  IXtalRecTool  reconstructs energy deposit intensity & centroid from digi information
 * -  IXtalDigiTool allows for the calculation of the digital response of one cal crystal given a list of enery depositions for that crystal. Includes trigger response.
 *
 * Currently there is only one concrete implementation of each of the CalResponseTools:  XtalRecTool, XtalDigiTool
 *
 * @subsection jobOptions jobOptions
 *
 * once again, suitable defaults exist for all parameters, all params are optional.
 * 
 * @param XtalEnergyTool.CalCalibSvc
 * where to retrieve calib data (defatul="CalCalibSvc")
 * @param XtalDigiTool.CalCalibSvc
 * where to retrieve calib data (defatul="CalCalibSvc")
 * @param XtalDigiTool.NoRandomNoise
 * disable calculation of random poissonic fluctuation & electronic noise
 * @param XtalDigiTool.SuperVerbose
 * enable extra debug logging.
 * 
 * 
 * @section CalXtalRecAlg CalXtalRecAlg
 *   CalXtalRecAlg takes the digitized calorimeter information from CalDigiCol
 *   as input, calculates the energy and position in each hitted crystal
 *   and stores this data into CalXtalRecCol.  CalXtalResponse package is 
 *   used for the estimation of energy & position from digi info.  See
 *   documentation in CalXtalResponse for details.
 * 
 * @subsection jobOptions jobOptions
 * @param CalXtalRecAlg.xtalRecTool
 *        name of CalXtalResponse/IXtalRecTool based tool performing xtal digi->energy conversion (default is "XtalRecTool")
 * @param CalXtalRecAlg.tupleFilename
  	  *        name of optional CalReconTuple output file. (default is "", which generates no tuple).
 *
 */

