// Mainpage for doxygen

/** @mainpage package CalUtil
 *
 * @authors R.Dubois, Zach Fewtrell
 *
 * @section description Description
 *
 * This package provides utilities for common use in the CAL.     
 *
 * @section CalFailureModeSvc CalFailureModeSvc
 * CalFailureModeSvc creates a list of large-scale failures in the CAL, 
 * and utilities
 * to search the lists to allow digi and recon algorithms to ignore hits based
 * on those lists.
 *
 * It can take lists of towers, (tower,AFEE) or (tower, Controller) pairs
 * to create the lists of dead objects. It provides a method to see if a 
 * given CalXtalId is contained in the lists.
 *
 * @subsection jobOptions jobOptions
 *
 * @param CalFailureModeSvc.towerList
 * Provide a list of towers that will be made dead.
 * Format: "tower"
 * @param CalFailureModeSvc.towerAfeeList
 * Provide a list of (tower, AFEE) pairs that will be made dead.
 * Format: "tower_afee". Afee runs 0-3: x+,y+,x-,y-
 * @param CalFailureModeSvc.towerControllerList
 * Provide a list of (tower, controller) pairs that will be made dead.
 * Format: "tower_controller". Controller runs 0-3: x+,y+,x-,y-.
 *
 *
 * @subsection CalDefs
 - CalDefs.h header file contains a wide array of contiguous index classes
   for indexing the various cal components (xtal, diode, xtal-face, adc-range, 
   etc.)
 - Internal integer representation is contiguous, so they are suitable for array
   indexing.
 - They have iteration operators defined so they are good for use in loops.  
 - They also have a full set of conversion routines defined for converting from
   one component type to another.  
 - These classes are wrappers for simple integers so they can be safely put on 
   function stack.  
 - Classes exist for indexing components within a single crystal or througou
   the entire LAT.  Appropriate conversion routines are provided.
 * 
 * @subsection CalVec
   - template wrapper for std::vector w/ the additional feature that it can only 
     be indexed by a single type.  Use of the index classes from CalDefs.h will
     allow for safe and convenient indexing of Cal based arrays.
   - commonly used std::vector functions are provided
 * 
 * @subsection CalArray
   - template class like CalVec, but uses internal array storage, instead of STL vector.
   - once again, has type checking on indexing to protect against accidental mis-indexing.
   - storage is allocated on construction for the proper size based on the index type selected.

 * @subsection unit_test
   - thoroughly tests the CalVec, CalArray & CalDefs classes.
   
 *
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr>
 * @section requirements requirements
 * @verbinclude requirements
 * <hr>
 *
 */

