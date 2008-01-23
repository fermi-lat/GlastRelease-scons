/** @mainpage package AcdUtil
* @author Heather Kelly
*
* @section intro Introduction
* This package contains all various ACD utility classes and functions.
* The common feature of the code in AcdUtil is that it is can be used
* anywhere else in the ACD software.  
*
* This tends to fall into two categories: 
*  1) Geomtrical descriptions and functions
*  2) Code to access and use calibrations
*
*
* This package contains three Gaudi services which are described more below
*  1) AcdGeometrySvc is used to get access to details of the AcdGeomerty
*  2) AcdCalibSvc provides the calibrations used in reconstrucion
*  3) AcdSimCalibSvc provides the calibrations used in simulation
*
*
* <hr>
* @section AcdGeometrySvc AcdGeometrySvc
*
* Provide access to detail of the acd geometry.  Mainly used to fill and give access
* to the AcdGeomMap, which contains all the AcdTileDim and AcdRibbonDim objects.
* Also provides a few utilities needed for getting the corners of the ACD and 
* keeping track of ACD elements. 
*
* @subsection jobOptions jobOptions
* None
* 
* @subsection AcdGeomMap AcdGeomMap
* The AcdGeomMap object gives information about all the ACD detector elements.
* Access is provided by idents::AcdId:
*  - const AcdRibbonDim* AcdGeomMap::GetRibbon(const idents::AcdId& id, IAcdGeometrySvc &detSvc)
*    - Get a geometrical description of an ACD ribbon
*  - const AcdTileDim* AcdGeomMap::GetTile(const idents::AcdId& id, IAcdGeometrySvc &detSvc)
*    - Get a geometrical description of an ACD tile
*
* AcdTileDim describes a tile as 1 or 2 flat planes, either rectangular or trapezoidal
*
* AcdRibbonDim describres a ribbons as a series of Rays.
*
* <hr>
* @section AcdCalibSvc AcdCalibSvc
*
* Provides access to  calibrations used in reconstruction.
* This service is used by the reconstruction code.  In particular, AcdRecon/AcdPha2MipTool.
* This service allows the user to select the flavor of each calibration seperately.
*
* There are two ways to set the calibration flavor.
* -# Setting the jobOption "DefaultFlavor" will change all the flavors
* -# Setting any of other Flavor* optiosn will change only that one flavor
*
* For safety, CalibDataSvc requires that all used flavors be registered.  Therefore, the user
* must insure that the flavor they select is also registers in CalibDataSvc jobOptions.
*
* @subsection jobOptions jobOptions
* @param CalibDataSvc
* The name of the service which provide the calibration data
*
* @param DefaultFlavor
* The flavor of calibration to use for any give calibration if not overridden
*
* @param FlavorPed
* The flavor of calibration to use for pedestals.
*
* @param FlavorGain
* The flavor of calibration to use for gains (aka MIP peaks).
*
* @param FlavorHighRange
* The flavor of calibration to use for the high range calibration.
*
* @param FlavorCoherentNoise
* The flavor of calibration to use for the coherent noise calibration.  
* 
*
*
* <hr>
* @section AcdSimCalibSvc AcdSimCalibSvc
*
* Provides access to  calibrations used in simulation.
* This service is used by the simulation code.  In particular, AcdDigi/AcdDigiAlg.
* This service allows the user to select the flavor of each calibration seperately.
*
* There are two ways to set the calibration flavor.
* -# Setting the jobOption "DefaultFlavor" will change all the flavors
* -# Setting any of other Flavor* optiosn will change only that one flavor
*
* For safety, CalibDataSvc requires that all used flavors be registered.  Therefore, the user
* must insure that the flavor they select is also registers in CalibDataSvc jobOptions.
*
* @subsection jobOptions jobOptions
* @param CalibDataSvc
* The name of the service which provide the calibration data
*
* @param DefaultFlavor
* The flavor of calibration to use for any give calibration if not overridden
*
* @param FlavorPed
* The flavor of calibration to use for pedestals.
*
* @param FlavorGain
* The flavor of calibration to use for gains (aka MIP peaks).
*
* @param FlavorVeto
* The flavor of calibration to use for veto thresholds.
*
* @param FlavorCno
* The flavor of calibration to use for CNO thresholds.
*
* @param FlavorRange
* The flavor of calibration to use for range crossover thresholds.
*
* @param FlavorHighRange
* The flavor of calibration to use for the high range calibration.
*
* @param FlavorCoherentNoise
* The flavor of calibration to use for the coherent noise calibration.  
*
*
* <hr>
* @section notes release notes
* release.notes
*
* @section requirements requirements
* @verbinclude requirements
*
*/

