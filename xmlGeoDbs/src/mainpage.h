// Mainpage for doxygen

/** @mainpage package xmlGeoDbs
 *
 * @author J. Bogart
 *
 * @section intro Introduction
 * This package provides a place to keep geometry description files
 * separately from the code (in packages xmlUtil and detModel) which
 * interprets them.  Files common to all descriptions (currently 
 * dtd and materials file) are kept in the xml directory of this package. 
 * Instrument-specific files are kept in a subdirectory of the xml 
 * directory named after the instrument so, for example, files specific
 * to the flight instrument are in xml/flight.  The top file is named
 * after the instrument (e.g., xml/flight/flight.xml).
 * The "meat" of the description is entirely in included files, such
 * as flightCALDimPrim.xml (floating point constants used in calorimeter
 * description) and flightTKROneTkr.xml (simple solids and their positionings
 * to build up tracker for a single tower).
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr>
 * @section requirements requirements
 * @verbinclude requirements
 * <hr> 
 * @todo  Reorganize so that top file flight.xml is just a shell
 * @todo  Define symbols in requirements for each instrument's top file
 */


