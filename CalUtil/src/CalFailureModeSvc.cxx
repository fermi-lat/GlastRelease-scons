// Implementation file for CalFailureModeSvc which handles the failure mode testing
// for the CAL.
// 
//
// $Header$
//
// Author: Richard Dubois


#include "CalUtil/CalFailureModeSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include <algorithm>


// declare the service factories for the CalFailureModeSvc
static SvcFactory<CalFailureModeSvc> a_factory;
const ISvcFactory& CalFailureModeSvcFactory = a_factory; 

CalFailureModeSvc::CalFailureModeSvc(const std::string& name,ISvcLocator* svc) : Service(name,svc)
{
    // Purpose and Method: Constructor - Declares and sets default properties
    //                     
    // Inputs: service name and locator 
    //         
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    // declare the properties

    declareProperty("towerList", m_towerListProperty);
    declareProperty("towerLayerList", m_towerLayerListProperty);
}

StatusCode  CalFailureModeSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ICalFailureModeSvc == riid) {
        *ppvIF = dynamic_cast<ICalFailureModeSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}

const IID&  CalFailureModeSvc::type () const {
    return IID_ICalFailureModeSvc;
}

StatusCode CalFailureModeSvc::initialize () 
{
    // Purpose and Method: Initialize the lists of dead units
    //                     
    // Inputs: None        
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    StatusCode  status = StatusCode::SUCCESS;
    
    // Open the message log
    MsgStream log( msgSvc(), name() );
    
    // Call super-class
    Service::initialize ();
    
    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }
    
    
    m_failureModes = 0;
    processTowerList();
    processTowerLayerList();
    
    return StatusCode::SUCCESS;
}

StatusCode CalFailureModeSvc::finalize () {return StatusCode::SUCCESS;}

void CalFailureModeSvc::processTowerLayerList() {
    // Purpose and Method: process the jobOptions input lists of (tower,layer) pairs
    //                     
    
    MsgStream log(msgSvc(), name());
    
    const std::vector<std::string>& theTowers = m_towerLayerListProperty.value( );
    if (theTowers.size() == 0) return;

    m_failureModes = m_failureModes || 1 << TOWER;

    log << MSG::DEBUG << "Towers and Layers to kill " << endreq;

    std::vector<std::string>::const_iterator it;
    std::vector<std::string>::const_iterator itend = theTowers.end( );
    for (it = theTowers.begin(); it != itend; it++) {
        int len = (*it).size();
        int delimPos = (*it).find_first_of('_');
        int tower = atoi((*it).substr(0, delimPos).c_str());
        int layer = atoi((*it).substr(delimPos+1, len-delimPos-1).c_str());
        log << MSG::DEBUG << "Tower " << tower << " Layer " << layer << endreq;
        std::vector<int>& curList = m_towerLayerList[tower];
        curList.push_back(layer);                
    }
}

void CalFailureModeSvc::processTowerList() {
     // Purpose and Method: process the jobOptions input lists of towers
    //                     
   
    MsgStream log(msgSvc(), name());
    
    const std::vector<std::string>& theTowers = m_towerListProperty.value( );

    if (theTowers.size() == 0) return;

    log << MSG::DEBUG << "Towers to kill " << endreq;

    m_failureModes = m_failureModes || 1 << TOWERLAYER;

    std::vector<std::string>::const_iterator it;
    std::vector<std::string>::const_iterator itend = theTowers.end( );
    for (it = theTowers.begin(); it != itend; it++) {
        int tower = atoi((*it).c_str());
        log << MSG::DEBUG << "Tower " << tower << endreq;
        m_towerList.push_back(tower);
    }
}

bool CalFailureModeSvc::matchTower(idents::CalXtalId id) {
    // Purpose and Method: check whether given id is in any of the identified lists
    //                     
    
    if (m_towerList.size() == 0) return false;

    int tower = id.getTower();
    
    // Search to see if this event id is among the list of ids we want to pause on
    int *loc = std::find(m_towerList.begin(), m_towerList.end(), tower);                

    return (loc != m_towerList.end());
}

bool CalFailureModeSvc::matchTowerLayer(idents::CalXtalId id) {
    // Purpose and Method: look for the given id in the tower list
    //                     
    
    if (m_towerLayerList.size() == 0) return false;

    int tower = id.getTower();
    int layer = id.getLayer();
    
    std::vector<int> &layerList = m_towerLayerList[tower];
    
    // Search to see if this (tower,layer) is among the list
    int *loc = std::find(layerList.begin(), layerList.end(), layer);                
    
    return (loc != layerList.end());
}

bool CalFailureModeSvc::matchChannel(idents::CalXtalId id) {
    // Purpose and Method: look for the given id in the (tower,layer) list
    //                     
     
    if (matchTower(id)) return true;
    if (matchTowerLayer(id)) return true;
    
    return false;
}
