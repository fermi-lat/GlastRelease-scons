
#ifndef _GlastSvc_IRandomAccess_H
#define _GlastSvc_IRandomAccess_H 1

// Include files
#include "GaudiKernel/IAlgTool.h"
#include <string>
#include <vector>

class HepRandomEngine;

/** @class IRandomAccess
* @brief Abstract definition of a tool to be called from the
* GlastRandomSvc to obtain the pointer to the local CLHEP static
* engine
* 
* $Header$
*/

static const InterfaceID IID_IRandomAccess("IRandomAccess", 1 , 0);

class IRandomAccess :virtual public IAlgTool {
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IRandomAccess;}

    virtual HepRandomEngine* getTheEngine()=0;
    virtual void setTheEngine(std::string engineName)=0;
};

#endif  // _GlastSvc_IRandomAccess_H
