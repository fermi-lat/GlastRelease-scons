// $Header$
// 
//  Original author: Toby Burnett tburnett@u.washington.edu
//

#include "FluxSvc/FluxSvc.h"

#include "Gaudi/Kernel/SvcFactory.h"
#include "Gaudi/MessageSvc/MsgStream.h"

#include "Gaudi/Kernel/Incident.h"
#include "Gaudi/Interfaces/IIncidentSvc.h"
#include "Gaudi/JobOptionsSvc/Property.h"

#include "Flux.h"

#include "flux/FluxMgr.h"
#include <algorithm>

// declare the service factories for the FluxSvc
static SvcFactory<FluxSvc> a_factory;
const ISvcFactory& FluxSvcFactory = a_factory;

void FATAL(const char* text){std::cerr << text << std::endl;}
void WARNING(const char* text){std::cerr << text << std::endl;}


// ------------------------------------------------
// Implementation of the FluxSvc class
// ------------------------------------------------
/// Standard Constructor
FluxSvc::FluxSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc)
{
    // declare the properties and set defaults
    
    //declareProperty ("size", m_default_name);

}

std::list<std::string> FluxSvc::fluxNames()const{
    return m_fluxMgr->sourceList();
}

StatusCode FluxSvc::source(std::string name, IFlux*& flux) {

    std::list<std::string> source_list( fluxNames() );

    if( std::find(source_list.begin(), source_list.end(), name) == source_list.end() )
        return StatusCode::FAILURE;
    
    flux =  new Flux(name);
    
    return StatusCode::SUCCESS;
}

/// Standard Destructor
FluxSvc::~FluxSvc()  
{
}


// initialize
StatusCode FluxSvc::initialize () 
{
    StatusCode  status =  Service::initialize ();
        
    // bind all of the properties for this service
    setProperties ();
    
    // open the message log
    MsgStream log( msgSvc(), name() );
 
    // create a FluxMgr object which will then be available.
    m_fluxMgr = new FluxMgr();

    Flux::mgr(m_fluxMgr); // tell our Flux object

    // check that it was made properly
    if( m_fluxMgr->sourceList().empty()) {
        log << MSG::ERROR  << "Did not initialize properly: no sources detected" << endreq;
        status = StatusCode::FAILURE;
    }


    return status;
}

// finalize
StatusCode FluxSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
       
    delete m_fluxMgr;
    return status;
}

