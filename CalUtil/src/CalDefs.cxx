// $Header $

/** @file
    @author fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/CalDefs.h"

// GLAST INCLUDES

// EXTLIB INCLUES

// STD INCLUDES
#include <string>
#include <sstream>

namespace CalUtil {
  using namespace std;

  const string DIR_MNEM[] = {
    "X",
    "Y"
  };

  const std::string &DirNum::toStr() const {
    return DIR_MNEM[m_data];
  }


  const string FACE_MNEM[] = {
    "POS",
    "NEG"
  };

  const std::string &FaceNum::toStr() const {
    return FACE_MNEM[m_data];
  }

  const string DIODE_MNEM[] = {
    "LRG",
    "SM"
  };

  const std::string &DiodeNum::toStr() const {
    return DIODE_MNEM[m_data];
  }


  const string THX_MNEM[] = {
    "X8",
    "X1"
  };

  const std::string &THXNum::toStr() const {
    return THX_MNEM[m_data];
  }


  const string RNG_MNEM[] = {
    "LEX8",
    "LEX1",
    "HEX8",
    "HEX1"
  };

  const std::string &RngNum::toStr() const {
    return RNG_MNEM[m_data];
  }


  const string ASYM_MNEM[] = {
    "ASYM_LL",
    "ASYM_LS",
    "ASYM_SL",
    "ASYM_SS"
  };


  /// generic type2string converter
  template <typename _T>
  string toString(const _T &val) {
    ostringstream tmp;
    tmp << val;
    return tmp.str();
  }

  const std::string &AsymType::toStr() const {
    return ASYM_MNEM[m_data];
  }


  std::string XtalIdx::toStr() const {
    return "T" + toString(getTwr().val()) +
      "L" + toString(getLyr().val()) +
      "C" + toString(getCol().val());
      
  }

  std::string FaceIdx::toStr() const {
    return getXtalIdx().toStr() + "F" + toString(getFace().val());
  }

  std::string DiodeIdx::toStr() const {
    return getFaceIdx().toStr() + "D" + toString(getDiode().val());
  }

  std::string RngIdx::toStr() const {
    return getFaceIdx().toStr() + "R" + toString(getRng().val());
  }

};
