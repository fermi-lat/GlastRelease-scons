#ifndef _IXtalPosTool_H
#define _IXtalPosTool_H 1
/*! @class IXtalPosTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for calculation of the deposited energy centroid for a single GLAST LAT calorimeter crystal.
 * 
 *
 */

// LOCAL INCLUDES

// GLAST INCLUDES
#include "idents/CalXtalId.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES

static const InterfaceID IID_XtalPosTool("IXtalPosTool", 1 , 0);

using namespace idents;

class IXtalPosTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_XtalPosTool; }

  /// calculate position of energy centroid in mm from xtal center, given the digital response for both faces 
  /// \param CalXtalId specify xtal log
  /// \param adcP input adc val for Positive end
  /// \param adcN input adc val for Negative end
  /// \param rngP input adc rng for Positive end
  /// \param rngN input adc rng for Negative end
  /// \param position output energy centroid position along xtal length in mm from center
  virtual StatusCode calculate(const CalXtalId &xtalId,
                               CalXtalId::AdcRange rngP,
                               CalXtalId::AdcRange rngN,
                               int adcP, 
                               int adcN, 
                               float &position
                               ) = 0;

};

#endif //_IXtalPosTool_H
