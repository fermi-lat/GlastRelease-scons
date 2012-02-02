
#ifndef HEPREPXMLSVC_H
#define HEPREPXMLSVC_H 1

/** 
 * @class HepRepXmlSvc
 *
 * @brief Load IRegister tools defined in HepRepSvc
 *
 * 
 * @author Heather Kelly 
 *
 * $Header$
 */

#include "HepRepSvc/IRegister.h"
#include "GaudiKernel/Service.h"

class HepRepXmlSvc : public Service
{
public:
    
    HepRepXmlSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~HepRepXmlSvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    
//    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

//    static const InterfaceID& interfaceID() {
//        return IHepRepXmlSvc::interfaceID(); 
//    }
    /// return the service type
//    const InterfaceID& type() const;

private:

    IRegister *m_xmlTool;


};

#endif // DETDISPLAYSVC_H
