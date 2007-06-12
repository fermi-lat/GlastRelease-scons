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
#include <vector>
#include <algorithm>

namespace {
  vector<string> tokenize_str(const string & str,
                              const string & delims)
  {
    // Skip delims at beginning.
    string::size_type lastPos = str.find_first_not_of(delims, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delims, lastPos);

    vector<string> tokens;

    while (string::npos != pos || string::npos != lastPos)
      {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delims.  Note the "not_of"
        lastPos = str.find_first_not_of(delims, pos);
        // Find next "non-delimiter"
        pos     = str.find_first_of(delims, lastPos);
      }

    return tokens;
  }

}

namespace CalUtil {
  using namespace std;

  const string DIR_MNEM[] = {
    "X",
    "Y"
  };

  const string &DirNum::toStr() const {
    return DIR_MNEM[m_data];
  }


  const string FACE_MNEM[] = {
    "POS",
    "NEG"
  };

  const string &FaceNum::toStr() const {
    return FACE_MNEM[m_data];
  }

  const string DIODE_MNEM[] = {
    "LRG",
    "SM"
  };

  const string &DiodeNum::toStr() const {
    return DIODE_MNEM[m_data];
  }

  ostream &operator<<(ostream &stream, const DiodeNum &id) {
      stream << id.toStr();
      
      return stream;
  }
  
  DiodeNum::DiodeNum(const std::string &str) {
    /// find string in list of mnemonics
    static const string *const first = DIODE_MNEM;
    static const string *const last = DIODE_MNEM + DiodeNum::N_VALS;
    const string *const pos = find(DIODE_MNEM,
                                   last,
                                   str);
    if (pos == last)
      throw runtime_error(str + " is inavlid DiodeNum string repr");
    
    /// set internal index value
    m_data = pos - first;
  }



  const string THX_MNEM[] = {
    "X8",
    "X1"
  };

  const string &THXNum::toStr() const {
    return THX_MNEM[m_data];
  }


  const string RNG_MNEM[] = {
    "LEX8",
    "LEX1",
    "HEX8",
    "HEX1"
  };

  const string &RngNum::toStr() const {
    return RNG_MNEM[m_data];
  }


  const string ASYM_MNEM[] = {
    "ASYM_LL",
    "ASYM_LS",
    "ASYM_SL",
    "ASYM_SS"
  };


  const string &AsymType::toStr() const {
    return ASYM_MNEM[m_data];
  }


  string XtalIdx::toStr() const {
    return "T" + toString(getTwr().val()) +
      "L" + toString(getLyr().val()) +
      "C" + toString(getCol().val());
      
  }

  XtalIdx::XtalIdx(const std::string &str) {
    /// extract individual fields from string
    vector<string> parts(tokenize_str(str,"TLC"));

    if (parts.size() != N_FIELDS)
      throw std::runtime_error("Invalid MeanDACZId string repr: " + str);
 
  }

  string FaceIdx::toStr() const {
    return getXtalIdx().toStr() + "F" + toString(getFace().val());
  }

  string DiodeIdx::toStr() const {
    return getFaceIdx().toStr() + "D" + toString(getDiode().val());
  }

  string RngIdx::toStr() const {
    return getFaceIdx().toStr() + "R" + toString(getRng().val());
  }
};
