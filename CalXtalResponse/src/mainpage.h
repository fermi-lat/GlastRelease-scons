// Mainpage for doxygen

/** @mainpage package CalXtalResponse

@section CalCalibSvc CalCalibSvc
 

ICalCalibSvc provides a Cal specific streamlined interface to the
GLAST calorimeter calibration database.  CalCalibSvc is currently 
the only concrete implementation of this interface.  

CalCalibSvc supports the following features:
 
- Simple interface: requires only specification of unique Cal xtal
and range as input.  

- Values are returned for the most part as C primitives (float, int,
etc..) C++ std::vectors are used where appropriate.

-  Gleam/Gaudi mechanics such as TDS access, validity checking,
calibration data format, database access, and data storage are
all transparent to the user.

- Some calibration types are vectors which represent the 'knots' on
a spline curve.  

CalCalibSvc generates the spline objects for these types for the
user.  CalCalibSvc handles storage/deallocation/validation for
these objects.  

Once created, spline objects are retained in memory until they
become invalid.  At this point, the cache is flushed and new
objects are created.

- CalCalibSvc can be configured to assign different flavors to each
calibration type.

- Multiple instances of CalCalibSvc may be created if the user
needs more than one flavor for one or more calibration types.


@subsection jobOptions jobOptions
 
Suitable default values are currently provided for all options, all
parameters are 'optional'.


@param CalibDataSvc
Name of service from which calibration data is
obtained. (default="CalibDataSvc")

@param DefaultFlavor
Flavor for all calibration types unless otherwise specified
(default="vanilla")

@param FlavorIntNonlin
calib flavor for integral nonlinearity constants.
@param FlavorAsym
override calib flavor for asymmetry data.
@param FlavorPed
override calib flavor for pedestal data.
@param FlavorMeVPerDac
override calib flavor for MeVPerDac data.
@param FlavorTholdCI
override calib flavor for tholdCI data.
@param FlavorTholdMuon
override calib flavor for tholdMuon data.

@section CalResponseTools CalResponseTools
 
The CalReponseTools are a collection of tools which provide a
consistent interface to calorimeter response calculations.  

CalResponseTools consists of 3 interfaces.
-  IXtalRecTool  reconstructs energy deposit intensity & centroid
from digi information

-  IXtalDigiTool allows for the calculation of the digital response
of one cal crystal given a list of enery depositions for that
crystal. Includes trigger response.

- ICalTrigTool generate FLE & FHE trigger response from digi info.
also populates GLtDigi TDS class & can run in either 1-range or 4-range
mode.  Can be called either per xtal or for entire Cal.
 
Currently there is only one concrete implementation of each of the
CalResponseTools:  XtalRecTool, XtalDigiTool, CalTrigTool

@subsection jobOptions jobOptions

once again, suitable defaults exist for all parameters, all params
are optional.

 
@param XtalEnergyTool.CalCalibSvc
where to retrieve calib data (defatul="CalCalibSvc")
@param XtalDigiTool.CalCalibSvc
where to retrieve calib data (defatult="CalCalibSvc")
@param CalTrigTool.CalCalibSvc
where to retrieve calib data (default="CalCalibSVc")
@param XtalDigiTool.NoRandomNoise
disable calculation of random poissonic fluctuation & electronic noise
 
@section CalXtalRecAlg CalXtalRecAlg
CalXtalRecAlg takes the digitized calorimeter information from CalDigiCol
as input, calculates the energy and position in each hitted crystal
and stores this data into CalXtalRecCol.  CalXtalResponse package is 
used for the estimation of energy & position from digi info.  See
documentation in CalXtalResponse for details.
 
@subsection jobOptions jobOptions
@param CalXtalRecAlg.xtalRecTool
name of CalXtalResponse/IXtalRecTool based tool performing
xtal digi->energy conversion (default is "XtalRecTool")

@param CalXtalRecAlg.tupleName
name of optional CalReconTuple 
file. (default is "", which generates no tuple).

@section unit_test
- CalXtalResponse unit_test thoroughly tests CalCalibSvc, XtalDigiTool
  XtalRecTool and CalTrigTool.
- Most tests are based on generation of Mc hits and validation of the 
  digitization and reconstruction of that input.
- By default the test runs on a subset of xtals, but it can also run
  on all crystals for a much longer, but much more comprehensive test 
  of a calibraiton set.

- Most tests are run w/ out simulated noise in digi in order to check the
precision of the code.  But special noise tests are run by rerunning the same 
input many times & checking that the distribution of the noise is correct.

- For each test, all didgital outputs are verified including LAT, FLE, FHE & ULD threholds.
Also, reconstructed energy and position are checked for accuracy.

- Samples from throughout the energy range of the instrument is checked and
  crystals are checked at several positions along their length.


*/

  
  
