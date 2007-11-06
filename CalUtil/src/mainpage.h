// Mainpage for doxygen

/** @mainpage package CalUtil
 *
 * @authors Zach Fewtrell
 *
 * @section description Description
 *
 * This package provides utilities for common use in the CAL.
 * Package intended for use by both Gaudi-based and non-Gaudi
 * based code
 *
 *
 * @subsection CalDefs
 - CalDefs.h header file contains a wide array of contiguous index classes
   for indexing the various cal components (xtal, diode, xtal-face, adc-range, 
   etc.) This is required because idents::CalXtalId uses non-contiguous bit-field structure.

   All CalDefs classes have the following properties:

 - Internal integer representation is contiguous, so they are suitable for array
   indexing.

 - Iteration operators defined for use in loops.  

 - Full set of conversion routines defined for converting between related indeces.
  
 - All classes are wrappers for simple integers, can efficiently be placed on function stack.

 - Classes exist for indexing components within a single crystal or througout
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

 * @subsection CalGeom
   - LAT Geometry related utilties

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

