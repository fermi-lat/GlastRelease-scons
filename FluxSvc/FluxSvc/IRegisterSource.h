#ifndef _H_IRegisterSource
#define _H_IRegisterSource

// includes
#include "GaudiKernel/IAlgTool.h"

class IFluxSvc;

static const InterfaceID IID_IRegisterSource("IRegisterSource", 1 , 0); 
 

/** @class IRegisterSource
 * @brief Abstract definition of a tool to be called from FluxSvc to load external ISpectrumFactory enteries
 * 
 * 
 * $Header$
*/

 class IRegisterSource : virtual public IAlgTool {
 public:
 /// Retrieve interface ID
 static const InterfaceID& interfaceID() { return IID_IRegisterSource; }
 // Actual operator function 
 virtual StatusCode registerMe(IFluxSvc*) = 0;
 };

#endif  // _H_IRegisterSource
