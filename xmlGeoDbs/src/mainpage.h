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
 * Subdirectories are used as follows:
 * <ul>
   <li><b>flight</b> is for standard LAT instrument.  There are two
        versions: flight.xml and flightSegVols.xml.  The only difference
        is that CsI crystals are built up as a stack of smaller CsI 
        volumes in the latter case </li>
   <li><b>latAssembly</b> has descriptions of the different missing-tower
          configurations as the LAT is being assembled.  Since there is
          only one configuration with a particular number of towers,
          the top-level files have names like <b>2TowerSegVols.xml</li>
   <li><b>em</b> is for the complete Engineering Model instrument,
       calorimeter + tracker </li>
   <li><b>minitower</b> is for the Engineering Model, tracker only </li>
   <li><b>cu</b> is for Calibration Unit</li>
   <li><b>Slab</b> trivial single-box geometry for testing</li>
   <li><b>flightSlab</b> integrates slab and current flight model</li>
   <li><b>misalignment</b>  mostly use standard flight geometry, but
       place towers one by one to produce known misalignments for
       testing. </li>
   <li><b>EM</b> Unused; obsolete. </li>
   </ul>
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr>
 * @section requirements requirements
 * @verbinclude requirements
 * <hr> 
 * @todo  Define symbols in requirements for each instrument's top file
 */


