// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include <cassert>

#include "FlightTest.h"
//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( FlightTest );

static const AlgFactory<FlightTest>  Factory;
const IAlgFactory& FlightTestFactory = Factory;

//------------------------------------------------------------------------
//! ctor
FlightTest::FlightTest(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
, m_detSvc(0), m_ini(0)
{
    // declare properties with setProperties calls
    declareProperty("tupleName",  m_tupleName="");
    
}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode FlightTest::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    if( m_tupleName.empty()) {log << MSG::ERROR << "tupleName property not set!"<<endreq;
        //return StatusCode::FAILURE;
    }
    // now try to find the GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc).isFailure()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    // get the ini file
    m_ini = const_cast<xml::IFile*>(m_detSvc->iniFile()); //OOPS!
    assert(4==m_ini->getInt("glast", "xNum"));  // simple check
    
   return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode FlightTest::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    /*! Causes the TDS to be searched, if the data is unavailable, the appropriate
    converter is called to retrieve the data from some persistent store, in this
    case from an IRF. 
    This call asks for the TdGlastData pointer to be given the variable name glastData.  
    We are asking the event data service, eventSvc(), to provide the data located at 
    "/Event/TdGlastData" which denotes the location of the glast detector data in 
    our Gaudi Event Model.
    */
    
    // Get the pointer to the event
    SmartDataPtr<TdGlastData> glastData(eventSvc(),"/Event/TdGlastData");
    
    // retrieve TKR data pointer for the event
    const SiData *tkrDigiData = glastData->getSiData();
    if (tkrDigiData == 0) {
        log << MSG::INFO << "No TKR Data available" << endreq;
    } else {
    }
    
    
    // retrieve CAL data pointer
    const CsIData *calDigiData = glastData->getCsIData();
    if (calDigiData == 0) {   
        log << MSG::INFO << "No CAL Data available" << endreq;
    } else {
    }
    
    // retrieve the ACD data pointer
    const IVetoData *acdDigiData = glastData->getVetoData();
    // Check to see if there is any ACD data for this event - if not provide a message
    if (acdDigiData == 0) 
        log << MSG::INFO << "No ACD Data available" << endreq;
    else {
    }
    
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FlightTest::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    return sc;
}

