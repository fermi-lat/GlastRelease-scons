#ifndef _ICalEnergyTool_H
#define _ICalEnergyTool_H 1
/*! @class ICalEnergyTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for calculation of total deposited energy for a single calorimeter crystal.
 * 
 *
 */

// // Include files
#include "GaudiKernel/IAlgTool.h"
#include "idents/CalXtalId.h"

static const InterfaceID IID_CalEnergyTool("CalEnergyTool", 1 , 0);

class ICalEnergyTool : virtual public IAlgTool {
public:
  static const InterfaceID& interfaceID() { return IID_CalEnergyTool; }

  /// estimate total deposited energy given the digital response for both faces
  /// \param xtalId specify xtal log
  /// \param adcP input adc value for Positive end
  /// \param adcN input adc value for Negative end
  /// \param rangeP input adc range for Positive end
  /// \param rangeN input adc range for Negative end
  /// \param position centriod of energy deposition
  /// \param energy output total deposited energy 
  virtual StatusCode calculate(const idents::CalXtalId &xtalId, 
                               idents::CalXtalId::AdcRange rangeP,
                               idents::CalXtalId::AdcRange rangeN,
                               int adcP, 
                               int adcN,
                               float position,
                               float &energy                    // output
                               ) = 0;
  
  /// estimate total deposited energy given the digital response for one face and a centroid-position 
  /// \param xtalId specify xtal log
  /// \param adc input adc value
  /// \param range input adc range
  /// \param position input energy centroid position
  /// \param energy output total deposited energy 
  virtual StatusCode calculate(const idents::CalXtalId &xtalId,
                               idents::CalXtalId::XtalFace face,
                               idents::CalXtalId::AdcRange range,
                               int adc, 
                               float position,
                               float &energy                    // output
                               ) = 0;


};

#endif //_ICalEnergyTool_H
