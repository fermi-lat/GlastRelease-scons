#ifndef _IXtalEneTool_H
#define _IXtalEneTool_H 1
/*! @class IXtalEneTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for calculation of total deposited energy for a single calorimeter crystal.
 * 
 *
 */

// LOCAL INCLUDES

// GLAST INCLUDES
#include "idents/CalXtalId.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES

using namespace idents;

static const InterfaceID IID_XtalEneTool("IXtalEneTool", 1 , 0);

class IXtalEneTool : virtual public IAlgTool {
 public:
  static const InterfaceID& interfaceID() { return IID_XtalEneTool; }

  /// estimate total deposited energy given the digital response for both faces
  /// \param CalXtalId specify xtal log
  /// \param adcP input adc value for Positive end
  /// \param adcN input adc value for Negative end
  /// \param rngP input adc rng for Positive end
  /// \param rngN input adc rng for Negative end
  /// \param position centriod of energy deposition
  /// \param energy output total deposited energy in MeV
  /// \param belowThresh returns TRUE if any LEX8 ADC is below 0.5 * LAC threshold
  virtual StatusCode calculate(const CalXtalId &xtalId, 
                               CalXtalId::AdcRange rngP,
                               CalXtalId::AdcRange rngN,
                               int adcP, 
                               int adcN,
                               float &energy,    //output
                               bool &belowThresh //output
                               ) = 0;

  /// estimate total deposited energy given the digital response for one face and a centroid-position 
  /// \param CalXtalId specify xtal log
  /// \param adc input adc value
  /// \param position input energy centroid position
  /// \param energy output total deposited energy 
  /// \param belowThresh returns TRUE if any LEX8 ADC is below 0.5 * LAC threshold
  virtual StatusCode calculate(const CalXtalId &xtalId,
                               int adc, 
                               float position,
                               float &energy,    // output
                               bool &belowThresh // output
                               ) = 0;

};

#endif //_IXtalEneTool_H
