
// Include files
#include "Gaudi/Kernel/SmartIF.h"
#include "Gaudi/Kernel/Bootstrap.h"
#include "Gaudi/Interfaces/IAppMgrUI.h"
#include "Gaudi/Interfaces/IProperty.h"
#include "Gaudi/JobOptionsSvc/Property.h"
#include "Gaudi/System/System.h"

extern "C" void DetDisplay_loadRef();
extern "C" void GuiSvc_loadRef();
//------------------------------------------------------------------------------
//
//  Package    : CalRecon
//
//  Description: Test Main Program
//
//------------------------------------------------------------------------------
#include <iostream>

int main( int argn, char** argc) {
    if(0) DetDisplay_loadRef();
    if(0) GuiSvc_loadRef();

    // get the path to this package from its root environment variable: if not there,
    // assume that we are in the root
    const char * local_path = ::getenv("DETDISPLAYROOT");
    std::string joboptions_file = std::string(local_path? local_path: "");

    joboptions_file +=  std::string("/src/test/jobOptions.txt");
    
    // Create an instance of an application manager
    IInterface* iface = Gaudi::createApplicationMgr();
    
    SmartIF<IProperty>     propMgr ( IID_IProperty, iface );
    SmartIF<IAppMgrUI>     appMgr  ( IID_IAppMgrUI, iface );
    
    // Set properties of algorithms and services
    StatusCode status = StatusCode::SUCCESS;
    if ( iface && propMgr == iface )    {
        status = propMgr->setProperty( StringProperty("JobOptionsPath", joboptions_file) );
    }
    else  {
        exit(-1);
    }

    
    // Run the application manager and process events
    if ( appMgr )   {
        status = appMgr->run();
    }
    else  {
        return 0;
    }
    
    // All done - exit
    return 0;
    
}
