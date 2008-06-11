#ifndef CalUtil_bit_util_h
#define CalUtil_bit_util_h

/** @file bit_util.h
    Generic bit manipulation utilities.
    @author fewtrell
 */

namespace CalUtil {
  /// specified bit (indexed from LSb=0) to specified value in
  /// specified integer reference
  ///
  /// @param datum reference to integer value to be modified
  /// INT_TYPE type of integer value to be modified
  /// @param bitNum index of bit to be modified (from lsb = 0)
  /// @param bitVal boolean value to set bit to
  template <typename INT_TYPE>
  void set_bit(INT_TYPE &datum,
               const unsigned short bitNum,
               const bool bitVal=true) {
    const INT_TYPE mask = 1 << bitNum;

    // CASE 1: bitVal = true
    if (bitVal) { 
      datum |= mask;
    } 

    // CASE 2: bitVal = false
    else { 
      datum &= ~mask;
    }
  }

  /// return boolean value of specified bit in specified integer
  template <typename INT_TYPE>
  bool check_bit(INT_TYPE &datum,
                 const unsigned short bitNum) {
    const INT_TYPE mask = 1 << bitNum;

    return datum & mask;
  }
}
#endif
