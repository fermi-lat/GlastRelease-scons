#ifndef _ICalPosTool_H
#define _ICalPosTool_H 1

// Include files
#include "GaudiKernel/IAlgTool.h"
#include "idents/CalXtalId.h"

static const InterfaceID IID_CalPosTool("CalReponsePos", 1 , 0);

class ICalPosTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_CalPosTool; }


  // calculate position given the digital response on both faces
  virtual StatusCode calculate(const idents::CalXtalId &xtalId,
                               int adcP, 
                               idents::CalXtalId::AdcRange rangeP,
                               int adcN, 
                               idents::CalXtalId::AdcRange rangeN,
                               float &position                // output
                               ) = 0;
  
};

#endif //_ICalPosTool_H
