#ifndef _ICalEnergyTool_H
#define _ICalEnergyTool_H 1

// // Include files
#include "GaudiKernel/IAlgTool.h"
#include "idents/CalXtalId.h"

static const InterfaceID IID_CalEnergyTool("CalEnergyTool", 1 , 0);

class ICalEnergyTool : virtual public IAlgTool {
public:
  static const InterfaceID& interfaceID() { return IID_CalEnergyTool; }

  virtual StatusCode calculate(const idents::CalXtalId &xtalId, 
                               idents::CalXtalId::AdcRange rangeP,
                               idents::CalXtalId::AdcRange rangeN,
                               int adcP, 
                               int adcN,
                               float position,
                               float &energy                    // output
                               ) = 0;
  
  // calculate energy from xtalId, one face/range/adc, and a position
  virtual StatusCode calculate(const idents::CalXtalId &xtalId,
                               idents::CalXtalId::AdcRange range,
                               idents::CalXtalId::XtalFace face,
                               int adc, 
                               float position,
                               float &energy                    // output
                               ) = 0;


};

#endif //_ICalEnergyTool_H
