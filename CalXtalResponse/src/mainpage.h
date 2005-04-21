// Mainpage for doxygen

/** @mainpage package CalXtalResponse
 * @section CalCalibSvc CalCalibSvc
 * 
 *
 * ICalCalibSvc provides a Cal specific streamlined interface to the GLAST calorimeter calibration
 * database.  CalCalibSvc is currently the only concrete implementation of this interface.  
 * CalCalibSvc supports the following features:
 * 
 * - Simple interface: requires only specification of unique Cal xtal and range as input.  
 * Values are returned for the most part as C primitives (float, int, etc..) C++ std::vectors are used where appropriate.
 * -  Gleam/Gaudi mechanics such as TDS access, validity checking, calibration data format, database access, and data storage are all transowner to the user.
 * - Some calibration types are vectors which represent the 'knots' on a spline curve.  
 * CalCalibSvc generates the spline objects for these types for the user.  CalCalibSvc handles storage/deallocation/validation for these objects.  
 * Once created, spline objects are retained in memory until their become invalid.  At this point, the cache is flushed and new objects are created.
 * - CalCalibSvc can be configured to assign different flavors to each calibration type.
 * - Multiple instances of CalCalibSvc may be created if the user needs more than one flavor for one or more calibration types.
 *
 * @subsection jobOptions jobOptions
 * 
 * Suitable default values are currently provided for all options, so consider tham all 'optional'.
 *
 * @param CalibDataSvc
 * Name of service from which calibration data is obtained. (default="CalibDataSvc")
 * @param DefaultFlavor
 * Flavor for all calibration types unless otherwise specified (default="vanilla")
 * @param FlavorGain
 * calib flavor for CAL_ElecGain constants
 * @param FlavorIntNonlin
 * calib flavor for integral nonlinearity constants.
 * @param FlavorLightAsym
 * calib flavor for light-asymetry data.
 * @param FlavorLightAtt
 * calib flavor for light-attenuation data.
 * @param FlavorMuSlope
 * calib flavor for muon-slope data.
 * @param FlavorPed
 * calib flavor for pedestal data.
 *
 * @section CalResponseTools CalResponseTools
 * 
 * The CalReponseTools are a collection of tools which provide a consistent interface to calorimeter response calculations.  
 * Once again, as much Gleam/Gaudi mechanics as possible are hidden from the user.
 *
 * CalResponseTools consists of 3 interfaces.
 * -  ICalEnergyTool allows for the calculation of deposited energy given cal-digi information for one cal crystal.
 * -  ICalPosTool allows for the calculation of the centroid position of the deposited energy given the cal-digi information for one cal crystal.
 * -  ICalAdcTool allows for the calculation of the digital response of one cal crystal given a list of enery depositions for that crystal.
 *
 * The <i>TestResponseTools</i> are the first implemenation of the CalResponseTools interface.  They are intended to duplicate as accurately as possible the calculations in CalRecon/CalDigi/CalUtil as of EngineeringModel release v3r0402p22.  They are intended for facilitating/testing the incorperation of the CalResponseTools interface in lieu of more sophisticated routines which should be forthcoming and which will use the same interface.
 *
 * -  TestEnergyTool implements ICalEnergyTool.  The code is lifted directly from CalXtalRecAlg in CalRecon v5r16p1.
 * -  TestPosTool implements ICalEnergyTool.  The code is lifted directly from CalXtalRecAlg in CalRecon v5r16p1.
 * -  TestAdcTool implelents ICalAdcTool.  The code is lifted from two places: CalDigiAlg in CalDigi/v1r3p7 and LinearConvertADC from CalUtil/v1r2p2.
 *
 *
 * @subsection jobOptions jobOptions
 *
 * once again, suitable defaults exist for all parameters, none are necesarily 'required'.
 * 
 * @param TestEnergyTool.startTime 
 * for use in selecting calib data with appropriate validity interval.
 * @param TestEnergyTool.calibFlavor
 * flavor of calib data to use.  "none" for simple 'dummy' constants.  (default="none")
 * @param TestPosTool.startTime 
 * for use in selecting calib data with appropriate validity interval.
 * @param TestPosTool.calibFlavor
 * flavor of calib data to use.  "none" for simple 'dummy' constants.  (default="none")
 * @param TestAdcTool.doFluctuations
 * approximate Poisson distribution by gaussian for numberElectrons. (default=true)
 * @param TestAdcTool.xmlFile
 * xml file for grabbing CalDigi configuration.  (default="$(CALDIGIROOT)/xml/CalDigi.xml")
 * 
 * 
 * @section CalXtalRecAlg CalXtalRecAlg
 *   CalXtalRecAlg takes the digitized calorimeter information from CalDigiCol
 *   as input, calculates the energy and position in each hitted crystal
 *   and stores this data into CalXtalRecCol.  CalXtalResponse package is 
 *   used for the estimation of energy & position from digi info.  See
 *   documentation in CalXtalResponse for details.
 *
 * 
 * @subsection jobOptions jobOptions
 * @param CalXtalRecAlg.xtalEneTool
 *        name of CalXtalResponse/IXtalEneTool based tool performing xtal digi->energy conversion (default is "XtalEneTool")
 * @param CalXtalRecAlg.xtalPosTool
 *        name of CalXtalResponse/IXtalPosTool based tool performing xtal digi->position conversion (default is "XtalPosTool")
 *
 */

