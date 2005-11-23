/** @file CurrentEventAlg.cxx
@brief CurrentEventAlg prints to an ASCII file the event and run id of the
       current event being processed

$Header$

*/

// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "facilities/Util.h"

#include <iostream>
#include <fstream>
#include <string>


/** 
* \class CurrentEventAlg
*
* \brief This is an Algorithm designed to set the event seeds by accessing a function member in GlastRandomSvc 
* \author Toby Burnett
* 
* $Header$
*/

class CurrentEventAlg : public Algorithm {
public:
    CurrentEventAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    
  std::string m_outputFileName;
  std::ofstream *m_eventFile;
};

static const AlgFactory<CurrentEventAlg>  Factory;
const IAlgFactory& CurrentEventAlgFactory = Factory;

CurrentEventAlg::CurrentEventAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator) 
{
    // declare properties with setProperties calls
    declareProperty("fileName", m_outputFileName="eventId.txt");
}

//! set parameters and attach to various perhaps useful services.
StatusCode CurrentEventAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    facilities::Util::expandEnvVar(&m_outputFileName);

    m_eventFile = new std::ofstream(m_outputFileName.c_str());
    if (!m_eventFile)  {
        MsgStream log( msgSvc(), name() );
        log << MSG::WARNING << "Failed to create output event file" << endreq;
        return StatusCode::FAILURE;
    }

   
    return sc;
}

//! process an event
StatusCode CurrentEventAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    
    if (!m_eventFile->good()) {
        MsgStream   log( msgSvc(), name() );
        log << MSG::WARNING << "Output file is bad, not writing" << endreq;
        return sc;
    }
    
    SmartDataPtr<Event::EventHeader> header(eventSvc(), EventModel::EventHeader);
    if(0==header) {
        MsgStream   log( msgSvc(), name() );
        log << MSG::WARNING << " could not find the event header" << endreq;
        return StatusCode::SUCCESS;
    }

    *m_eventFile << "run=" << header->run() << " event=" << header->event() << "  ";

    m_eventFile->flush();
    return sc; 
}

StatusCode CurrentEventAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    m_eventFile->close();
    delete m_eventFile;
    return sc;
}

