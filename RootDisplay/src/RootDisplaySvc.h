
#ifndef ROOTDISPLAYSVC_H
#define ROOTDISPLAYSVC_H 1

/** 
 * @class RootDisplaySvc
 *
 * @brief Load IGuiTools 
 *
 * 
 * @author Heather Kelly 
 *
 * $Header$
 */

#include "GuiSvc/IGuiTool.h"
#include "GaudiKernel/Service.h"

class RootDisplaySvc : public Service
{
public:
    
    RootDisplaySvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~RootDisplaySvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    
    //StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

//    static const InterfaceID& interfaceID() {
//        return IRootDisplaySvc::interfaceID(); 
//    }
    /// return the service type
//    const InterfaceID& type() const;

private:

    IGuiTool *m_rootMenuTool;


};

#endif // ROOTDISPLAYSVC_H
