#ifndef _ICalAdcTool_H
#define _ICalAdcTool_H 1

// Include files
#include "GaudiKernel/IAlgTool.h"

#include "idents/CalXtalId.h"
#include "Event/MonteCarlo/McIntegratingHit.h"

#include <vector>

static const InterfaceID IID_CalAdcTool("CalReponseAdc", 1 , 0);

class ICalAdcTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_CalAdcTool; }

  virtual StatusCode calculate(const idents::CalXtalId &xtalId, 
                               const std::vector<const Event::McIntegratingHit*> &hitList,
                               std::vector<int> &adcP,              // output - ADC's for all ranges 0-3
                               idents::CalXtalId::AdcRange &rangeP, // output - best range
                               std::vector<int> &adcN,              // output - ADC's for all ranges 0-3
                               idents::CalXtalId::AdcRange &rangeN  // output - best range
                               ) = 0;


};

#endif //_ICalAdcTool_H
