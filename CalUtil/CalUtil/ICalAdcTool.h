#ifndef _ICalAdcTool_H
#define _ICalAdcTool_H 1
/*! @class ICalAdcTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for estimation of the digital response of one calorimeter crystal to a set of energy depositions.
 * 
 *
 */


// Include files
#include "GaudiKernel/IAlgTool.h"

#include "idents/CalXtalId.h"
#include "Event/MonteCarlo/McIntegratingHit.h"

#include <vector>

static const InterfaceID IID_CalAdcTool("CalReponseAdc", 1 , 0);

class ICalAdcTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_CalAdcTool; }

  /// calculate Adc response for one cal xtal.  also select best range.
  /// \param xtalId specify xtal log
  /// \param hitList input vector of energy depositions.  (const *) is used to save space.
  /// \param rangeP output best range for Positive xtal face
  /// \param rangeN output best range for Negative xtal face
  /// \param adcP output vector of 4 adc responses for each range for Positive face - ranges always in default order (from 0-3)
  /// \param adcN output vector of 4 adc responses for each range for Negative face - ranges always in default order (from 0-3)
  /// \param lacP output boolean for log accept on Positive xtal face
  /// \param lacN output boolean for log accept on Negative xtal face
  virtual StatusCode calculate(const idents::CalXtalId &xtalId, 
                               const std::vector<const Event::McIntegratingHit*> &hitList,
                               bool &lacP,                          // output - log accept.
                               bool &lacN,                          // output - log accept. 
                               idents::CalXtalId::AdcRange &rangeP, // output - best range
                               idents::CalXtalId::AdcRange &rangeN, // output - best range
                               std::vector<int> &adcP,              // output - ADC's for all ranges 0-3
                               std::vector<int> &adcN               // output - ADC's for all ranges 0-3
                               ) = 0;


};

#endif //_ICalAdcTool_H
