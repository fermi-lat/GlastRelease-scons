#ifndef IXtalDigiTool_H
#define IXtalDigiTool_H
// $Header$
/** @file
    @author Z.Fewtrell
*/
/*! @class IXtalDigiTool
 *
 * \author Z.Fewtrell
 *
 * \brief Abstract interface class for estimation of the digital response of 
 one calorimeter crystal to a set of energy depositions.
 * 
 *
 */

// LOCAL INCLUDES
#include "ICalSignalTool.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES

namespace Event {
  class CalDigi;
}

static const InterfaceID IID_IXtalDigiTool("IXtalDigiTool", 1, 3);

class IXtalDigiTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_IXtalDigiTool; }

  /// represent all 4 diode signals in single xtal.
  typedef CalUtil::CalVec<CalUtil::XtalDiode, float> XtalSignalMap;

  /** \brief calculate Adc response for one cal xtalIdx.  also select best rng.

  \param cidac   input electronic signal level for each diode in crystal
  \param calDigi output empty CalDigi object to be populated (xtalId & rangeMode expected populated)
  \param lacBits output boolean for log accept on each xtal face
  \param zeroSuppress.  if zero suppression is on, i can optimize by not fully evaluating xtals which will not be recorded.

  */

  virtual StatusCode calculate(Event::CalDigi &calDigi,
                               CalUtil::CalVec<CalUtil::FaceNum, bool> &lacBits,
                               bool zeroSuppress, string calFirstRng= "autoRng"
                               ) = 0;

};

#endif //_IXtalDigiTool_H
