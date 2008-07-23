/** @file
    Mainpage for doxygen
    @author Z.Fewtrell
*/
// $Header$

/** 
    @mainpage package CalXtalResponse
    
    @section CalCalibSvc CalCalibSvc
    
    
    ICalCalibSvc provides a Cal specific streamlined interface to the
    GLAST calorimeter calibration database.  CalCalibSvc is currently 
    the only concrete implementation of this interface.  
    
    CalCalibSvc supports the following features:
    
    - Simple interface: requires only specification of unique Cal xtal
    and range as input.  
    
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
    
    @section CalResponseTools CalResponseTools
    
    The CalReponseTools are a collection of tools which provide a
    consistent interface to calorimeter response calculations.  
    
    CalResponseTools consists of 5 interfaces.
    - ICalSignalTool - calculate diode signal response (in charge injection dac (CIDAC) units) 
    for each crystal in Cal using TDS McIntegratingHitCol and IXtalSignalTool

    - IXtalSignalTool - calculate signal response for single crystal / McIntegrating hit pair.

    -  IXtalDigiTool allows for the calculation of the digital response
    of one cal crystal from diode signal levels.

    - ICalTrigTool generate FLE & FHE trigger response from digi or MC
    TDS data.

    -  IXtalRecTool  reconstructs energy deposit intensity & centroid
    from digi information
    
    
    
    Currently there is only one concrete implementation of each of the
    CalResponseTools:  CalSignalTool, XtalSignalTool, XtalDigiTool, CalTrigTool, XtalRecTool
    
    @subsection jobOptions jobOptions
    
    once again, suitable defaults exist for all parameters, all params
    are optional. see individual tools doxygen page for jobOptions
    
    @section CalXtalRecAlg CalXtalRecAlg
    CalXtalRecAlg takes the digitized calorimeter information from CalDigiCol
    as input, calculates the energy and position in each hit crystal
    and stores this data into CalXtalRecCol.  CalXtalResponse package is 
    used for the estimation of energy & position from digi info.  See
    documentation in CalXtalResponse for details.
    
    @subsection jobOptions jobOptions
    @param CalXtalRecAlg.XtalRecToolName
    name of CalXtalResponse/IXtalRecTool based tool performing
    xtal digi->energy conversion (default is "XtalRecTool")
    
    
    @section CalTupleAlg CalTupleAlg
    
    CalTupleAlg generates CalTuple entries from TDS digi data.  Current entries
    include:
    - CalXtalAdcPed: pedestal subtraced adc values per channel. (BESTRANGE)
    - CalXtalAdcRng: adc range selection per channel
    - CalXtalFaceSignal: Signal @ each crystal face in units of MeV deposited at
    center of xtal

    Each of these branches is a multi-dimensional array matching 
    CAL geometry ([16][8][12][2],  representing tower, layer, 
    xtal & xtal face respectively)

    - CalXtalAdcPedAllRange: pedestal subtracted adc values 
    for all adc channels (may be zero if data is not available).
    shape is [16][8][12][2][4] (last index represents ADC range.
    - CalXtalFaceSignalAllRange: face signal
    for all adc channels (may be zero if data is not available).
    shape is [16][8][12][2][4] (last index represents ADC range.
    
    @subsection jobOptions jobOptions
    
    @param CalTupleAlg.tupleName override name of CalTuple tree. 
    (default is "CalTuple")
    
    @param CalTupleAlg.tupleFilename optional name of CalTuple file.
    instructs ntupleWriterSvc to create CalTuple in it's own file
    instead of sharing the default file w/ other modules.  The
    default ("") will use the shared file
    
    @param CalCalibSvc CalCalibSvc
    specifies which ICalCalibSvc object should be used by CalTupleAlg.  
    default is "CalCalibSvc"
    
    @section CalFailureModeSvc CalFailureModeSvc
    CalFailureModeSvc creates a list of large-scale failures in the CAL, 
    and utilities
    to search the lists to allow digi and recon algorithms to ignore hits based
    on those lists.
    
    It can take lists of towers, (tower,AFEE) or (tower, Controller) pairs
    to create the lists of dead objects. It provides a method to see if a 
    given CalXtalId is contained in the lists.
    
    @subsection jobOptions jobOptions
    
    @param CalFailureModeSvc.towerList
    Provide a list of towers that will be made dead.
    Format: "tower"
    @param CalFailureModeSvc.towerAfeeList
    Provide a list of (tower, AFEE) pairs that will be made dead.
    Format: "tower_afee". Afee runs 0-3: x+,y+,x-,y-
    @param CalFailureModeSvc.towerControllerList
    Provide a list of (tower, controller) pairs that will be made dead.
    Format: "tower_controller". Controller runs 0-3: x+,y+,x-,y-.
    
    
    @section unit_test
    - CalXtalResponse unit_test thoroughly tests CalCalibSvc, XtalDigiTool
    XtalRecTool and CalTrigTool.
    - Most tests are based on generation of Mc hits and validation of the 
    digitization and reconstruction of that input.
    - By default the test runs on a subset of xtals, but it can also run
    on all crystals for a much longer, but much more comprehensive test 
    of a calibraiton set.
    
    - Most tests are run w/ out simulated noise in digi in order to check the
    precision of the code.  But special noise tests are run by rerunning the 
    same input many times & checking that the distribution of the noise is 
    correct.
    
    - For each test, all didgital outputs are verified including LAT, FLE, FHE &
    ULD threholds.  Also, reconstructed energy and position are checked for 
    accuracy.
    
    - Samples from throughout the energy range of the instrument is checked and
    crystals are checked at several positions along their length.
    
*/
