/** @file implement CalDiagnosticWord
    @author Z. Fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/CalDiagnosticWord.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/bit_util.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES

namespace CalUtil {
  /// build integer bitfields from local boolean arrays
  /// in accordance with "TEM Programming ICD Specification", LATDoc ID: LAT-XR-07314-01
  unsigned CalDiagnosticWord::getDatum() const {
    unsigned datum = 0;

    for (FaceNum face; face.isValid(); face++) {

      /// 32 bit word is divided between POS  & NEG face.
      /// first, bit of data for each face
      const unsigned short faceFirstBit = (face == POS_FACE) ? 0 : 16;

      // LAC BITS (fill bits 0-11 of 16 bit face word)
      for (ColNum col; col.isValid(); col++) {
        const unsigned short bitNum = faceFirstBit + col.val();
        set_bit(datum, bitNum, m_lacBits[face][col]);
      }

      
      // TRIG BITS (fill bits 12 & 13 of 16 bit face word)
      for (DiodeNum diode; diode.isValid(); diode++) {
        // 12 + 0 for LE, 1 for HE
        const unsigned short bitNum = faceFirstBit + ColNum::N_VALS + diode.val();

        set_bit(datum, bitNum, m_trigBits[face][diode]);
      }

    }

    return datum;
  }

  CalDiagnosticWord::CalDiagnosticWord(const unsigned dataWord) {
    /// populate my boolean arrays from bits in dataWord
    for (FaceNum face; face.isValid(); face++) {
      /// 32 bit word is divided between POS  & NEG face.
      /// first, bit of data for each face
      const unsigned short faceFirstBit = (face == POS_FACE) ? 0 : 16;

      // LAC BITS (fill bits 0-11 of 16 bit face word)
      for (ColNum col; col.isValid(); col++) {
        const unsigned short bitNum = faceFirstBit + col.val();
        m_lacBits[face][col] = check_bit(dataWord, bitNum);
      }

      // TRIG BITS (fill bits 12 & 13 of 16 bit face word)
      for (DiodeNum diode; diode.isValid(); diode++) {
        // 12 + 0 for LE, 1 for HE
        const unsigned short bitNum = faceFirstBit + ColNum::N_VALS + diode.val();
        
         m_trigBits[face][diode] = check_bit(dataWord, bitNum);
      }
    }
  }
};
