// Implementation file for AcdFailureModeSvc which handles the failure mode testing
// for the ACD.
// 
//
// $Header$
//
/** @file     
@author H.Kelly
*/


#include "AcdUtil/AcdFailureModeSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include <algorithm>


// declare the service factories for the AcdFailureModeSvc

static SvcFactory<AcdFailureModeSvc> a_factory;
const ISvcFactory& AcdFailureModeSvcFactory = a_factory; 

AcdFailureModeSvc::AcdFailureModeSvc(const std::string& name,ISvcLocator* svc) : Service(name,svc)
{
    // Purpose and Method: Constructor - Declares and sets default properties
    //                     
    // Inputs: service name and locator 
    //         
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None


    // declare the properties        
    declareProperty("detectorList", m_detectorListProperty);
}

StatusCode  AcdFailureModeSvc::queryInterface (const InterfaceID& riid, void **ppvIF) {

    if (IID_IAcdFailureModeSvc == riid) {
        *ppvIF = dynamic_cast<IAcdFailureModeSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}


StatusCode AcdFailureModeSvc::initialize () 
{
    // Purpose and Method: Initialize the lists of dead units
    //                     
    // Inputs: None        
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None


    StatusCode  status = StatusCode::SUCCESS;

    // Open the message log
    MsgStream msglog(msgSvc(), name());   
    msglog << MSG::INFO << "initialize" << endreq;

    // Call super-class
    Service::initialize ();

    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        msglog << MSG::ERROR << "Failed to set properties" << endreq;
    }

    processDetectorList();

    return StatusCode::SUCCESS;
}


StatusCode AcdFailureModeSvc::finalize () {return StatusCode::SUCCESS;}


void AcdFailureModeSvc::processDetectorList() {

    // Purpose and Method: process the jobOptions input lists of tiles
    //                     

    MsgStream msglog(msgSvc(), name());

    const std::vector<unsigned int>& theDetectors = m_detectorListProperty.value( );
    msglog << MSG::DEBUG << "Number of Detectors to kill " 
           << theDetectors.size() << endreq;

    if (theDetectors.size() == 0) return;

    msglog << MSG::INFO << "ACD Detectors to kill " << endreq;

    std::vector<unsigned int>::const_iterator it;
    std::vector<unsigned int>::const_iterator itend = theDetectors.end( );

    for (it = theDetectors.begin(); it != itend; it++) {
        //unsigned int detector = atoi((*it).c_str());
        unsigned int detector = *it;

        msglog << MSG::INFO << "Detector " << detector << endreq;

        m_detectorList.push_back(detector);
    }
}


bool AcdFailureModeSvc::matchAcdId(idents::AcdId id) {

    // Purpose and Method: check whether given id is in the list of detectors 

    if (m_detectorList.size() == 0) return false;

    unsigned int detector = id.id();

    // Search to see if this event id is among the list of ids we want to pause on

    std::vector<unsigned int>::iterator loc = std::find(m_detectorList.begin(), 
        m_detectorList.end(), detector);
    return (loc != m_detectorList.end());
}



