// $Header$

// Include files
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/IProperty.h"
#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/System.h"
//------------------------------------------------------------------------------
//
//  Package    : GuiSvc
//
//  Description: Graphics Main Program
//
//------------------------------------------------------------------------------
/**
Standard Glast-Gaudi main.
Specify the job options file location by:
1) specification, local or global, on command line.
2) content of env var JOBOPTIONS
3) src/jobOptions

*/
// declare function that may reduce the execution priority
void setPriority();

int main( int argn, char** argc) {
    
    std::string joboptions_file="src/jobOptions.txt"; // default
    
    const char* job = ::getenv("JOBOPTIONS"); // check for env var
    
    if( argn>1 ) { joboptions_file = argc[1];} // priority to command arg.
    else if( job ) { joboptions_file = job; }
    std::cerr << "Starting Glast-Gaudi job with job options file " 
        << joboptions_file << std::endl;
    
    
    // Create an instance of an application manager
    IInterface* iface = Gaudi::createApplicationMgr();
    
    SmartIF<IProperty>     propMgr ( IID_IProperty, iface );
    SmartIF<IAppMgrUI>     appMgr  ( IID_IAppMgrUI, iface );
    
    // Set properties of algorithms and services
    StatusCode status = StatusCode::SUCCESS;
    if ( iface && propMgr == iface )    {
        status = propMgr->setProperty( "JobOptionsPath", joboptions_file );
    }
    else  {
        exit(-1);
    }
    
    // Run the application manager and process events
    if ( appMgr )   {
        setPriority();
        status = appMgr->run();
    }
    else  {
        return 0;
    }
    // All done - exit
    return 0;
    
}
