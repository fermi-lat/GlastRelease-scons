#ifndef _H_IGuiSvc_
#define _H_IGuiSvc_

// includes
#include "GaudiKernel/IInterface.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IGuiSvc(900, 1 , 1); 

//forward declarations
namespace gui { class GuiMgr; }

/*! Proper Gaudi abstract interface class for the GuiSvc 
*/
class IGuiSvc : virtual public IInterface
{  

public:
    virtual gui::GuiMgr* guiMgr()=0;

    /// set flag to pause a run
    virtual void pause()=0;

    virtual void quit()=0;
    virtual void finishSetup()=0;  // dummy

    virtual StatusCode initialize ()=0;
    
    virtual StatusCode finalize ()=0;


    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IGuiSvc; }

};


#endif // _H_IGuiSvc
