#ifndef CalCalibDefs_h
#define CalCalibDefs_h

// LOCAL

// GLAST
#include "CalUtil/CalDefs.h"

// EXTLIB

// STD

/** @file 
    @author Zach Fewtrell
    
    @brief Shared definitions used by CalCalibSvc
    
*/

using namespace std;
using namespace CalUtil;
using namespace idents;

namespace CalXtalResponse {
  /** \brief id for 4 classes of Cal xtal asymmetry 

  4 asymmetry classes are
  - LL = large diode on both faces
  - LS = large diode on positive face, small diode on neg
  - SL = sm. diode on pos. face, large diode on neg
  - SS = small diode on both faces
      
  */

  class AsymType : public SimpleId {
  public:
    AsymType() : SimpleId() {}
    AsymType(DiodeNum posDiode, DiodeNum negDiode) :
      SimpleId(posDiode.val()*2 + negDiode.val()) {};
    
    inline DiodeNum getDiode(FaceNum face) {
      switch ((CalXtalId::XtalFace)face) {
      case CalXtalId::POS:
        return DiodeNum(m_data/2);
      case CalXtalId::NEG:
        return DiodeNum(m_data%2);
      default:
        throw invalid_argument("Bad FaceNum value. Programmer error");
      }
    }

    static const unsigned short N_VALS = 4; 

    bool isValid() const {return m_data < N_VALS;}


    bool operator==(const AsymType &that) const {return m_data == that.m_data;}
    bool operator!=(const AsymType &that) const {return m_data != that.m_data;}

  };
  const AsymType ASYM_LL(LRG_DIODE, LRG_DIODE);
  const AsymType ASYM_LS(LRG_DIODE, SM_DIODE);
  const AsymType ASYM_SL(SM_DIODE, LRG_DIODE);
  const AsymType ASYM_SS(SM_DIODE, SM_DIODE);
};

#endif
