// Mainpage for doxygen

/** @mainpage package calibUtil
 *
 * @authors J. Bogart, M. Gulati
 *
 * @section intro Introduction
 * This package contains (well, will contain) utilities for storing,
 * searching for, and accessing calibration data.  "Calibration data"
 * include at a minimum hardware status information, such as hot and
 * dead tracker strips, and conversion constants, such as those 
 * characterizing light attenuation in calorimeter crystals.
 *
 * The overall strategy is to store a limited amount of <b>metadata</b> 
 * about each set of calibration data in a relational database table, 
 * while the data are kept in flat files, one per result set.  Most 
 * likely two formats will be supported for these files: XML and ROOT.
 * See the <a href="http://www.slac.stanford.edu/exp/glast/ground/software/calibration/CalibSvcSpec.shtml">Calibration Services Specification</a> document
 * for more details.
 * 
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr>
 * @section requirements requirements
 * @verbinclude requirements
 * <hr> 
 * @todo  practically everything.
 *
 */

