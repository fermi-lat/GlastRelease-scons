#ifndef CALIBDATASVC_IINSTRUMENTNAME_H
#define CALIBDATASVC_IINSTRUMENTNAME_H

#include "GaudiKernel/IInterface.h"

/**  @class IInstrumentName
     @brief Simple interface to keep track of which instrument (LAT,
     EM, etc.) the process is concerned with.  Modeled after
     IDetDataSvc handling of event time.

     Intention is to implement in CalibDataSvc.
*/
static const InterfaceID IID_IInstrumentName(??, 1, 0);
class IInstrumentName {
public:

  virtual const bool validInstrumentName() = 0;
  virtual const std::string& getInstrumentName() = 0;
  virtual void setInstrumentName(const std::string& name) = 0;
}

#endif
