#ifndef _ICalPosTool_H
#define _ICalPosTool_H 1
/*! @class ICalPosTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for calculation of the deposited energy centroid for a single GLAST LAT calorimeter crystal.
 * 
 *
 */

// Include files
#include "GaudiKernel/IAlgTool.h"
#include "idents/CalXtalId.h"

static const InterfaceID IID_CalPosTool("CalReponsePos", 1 , 0);

class ICalPosTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_CalPosTool; }


  /// calculate position of energy centroid relative to xtal, given the digital response for both faces 
  /// \param xtalId specify xtal log
  /// \param adcP input adc value for Positive end
  /// \param adcN input adc value for Negative end
  /// \param rangeP input adc range for Positive end
  /// \param rangeN input adc range for Negative end
  /// \param position output energy centroid position
  virtual StatusCode calculate(const idents::CalXtalId &xtalId,
                               int adcP, 
                               idents::CalXtalId::AdcRange rangeP,
                               int adcN, 
                               idents::CalXtalId::AdcRange rangeN,
                               float &position
                               ) = 0;
  
};

#endif //_ICalPosTool_H
