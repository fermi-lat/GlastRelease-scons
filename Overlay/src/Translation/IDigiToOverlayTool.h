/** @file IDigiToOverlayTool.h

    @brief declaration of the IDigiToOverlayTool class

$Header$

*/

#ifndef IDigiToOverlayTool_h
#define IDigiToOverlayTool_h

#include "GaudiKernel/IAlgTool.h"

/** @class IBackgroundBintTool
    @brief Interface to tools to determine bin for overlay input
*/

static const InterfaceID IID_IDigiToOverlayTool("IDigiToOverlayTool", 1 , 0);

class IDigiToOverlayTool : virtual public IAlgTool
{
public:

    // Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IDigiToOverlayTool; }

    ///! The current value of the quantity that we are selecting on
    virtual StatusCode translate() = 0;
};


#endif
