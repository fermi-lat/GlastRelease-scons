#ifndef _H_IGuiSvc_
#define _H_IGuiSvc_

// includes
#include "Gaudi/Interfaces/IInterface.h"

//forward declarations
namespace gui { class GuiMgr; }

/*! Proper Gaudi abstract interface class for the GuiSvc 
*/
class IGuiSvc : virtual public IInterface
{  

public:
    virtual gui::GuiMgr* guiMgr()=0;

    virtual void quit()=0;
    virtual void finishSetup()=0;  // dummy

    virtual StatusCode initialize ()=0;
    
    virtual StatusCode finalize ()=0;
};


#endif // _H_IGuiSvc
