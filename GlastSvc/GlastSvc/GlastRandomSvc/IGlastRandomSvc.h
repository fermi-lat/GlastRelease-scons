
#ifndef _GlastSvc_IGlastRandomSvc_H
#define _GlastSvc_IGlastRandomSvc_H 1

#include "GaudiKernel/IInterface.h"
#include <vector>
#include <string>

// Declaration of the interface ID ( interface id, major version,
// minor version)

static const InterfaceID IID_IGlastRandomSvc("IGlastRandomSvc", 1 , 0);

// static const InterfaceID IID_IRandomAccess(912 , 1 , 0);

/** @class IGlastRandomSvc
* @brief Abstract definition of a service to set the CLHEP
*       random engine and seeds for Dll's that need random numbers
* @author Karl Young
*
* $Header:
*/

class IGlastRandomSvc : virtual public IInterface {
public:    

  /// initialization includes getting random engine addresses
  /// and setting engine type to use 
  virtual StatusCode initialize ()=0;

  /// perform cleanup  
  virtual StatusCode finalize ()=0;

  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_IGlastRandomSvc; }

};

#endif // _GlastSvc_IGlastRandomSvc_H




