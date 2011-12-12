// $Header$

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Overlay/IOverlayDataSvc.h"

#include "OverlayEvent/OverlayEventModel.h"
#include "OverlayEvent/EventOverlay.h"
#include "OverlayEvent/AcdOverlay.h"
#include "OverlayEvent/CalOverlay.h"
#include "OverlayEvent/TkrOverlay.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include <string>

/*! \class OverlayTestAlg
\brief 
  This is a place to put test code to examine what G4Generator has done.
  */

class OverlayTestAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    OverlayTestAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    IOverlayDataSvc* m_overlayInputSvc;
};


//static const AlgFactory<OverlayTestAlg>  Factory;
//const IAlgFactory& OverlayTestAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(OverlayTestAlg);

//
OverlayTestAlg::OverlayTestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)
{
//    declareProperty("StartX", m_startx = 200.); 
//    declareProperty("StartX", m_starty = 200.); 
}


/*! */
StatusCode OverlayTestAlg::initialize() 
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initializing..." << endreq;
    StatusCode  sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // Try to retrieve the input data service
    IService* tmpService = 0;
    if (service("OverlayInputSvc", tmpService, true).isFailure())
    {
        log << MSG::INFO << "No OverlayInputSvc available, no input conversion will be performed" << endreq;
    }
    else 
    {
        m_overlayInputSvc = SmartIF<IOverlayDataSvc>(tmpService);
        //m_overlayInputSvc = SmartIF<IOverlayDataSvc>(IID_IOverlayDataSvc, tmpService);

        if (!m_overlayInputSvc) log << MSG::INFO << "Could not find input data service interface" << endreq;
        else                  log << MSG::INFO << "Input data service successfully retrieved" << endreq;
    }

    // Now look for the output data service
    if (service("OverlayOutputSvc", tmpService).isFailure())
    {
        log << MSG::INFO << "No OverlayOutputSvc available, no output conversion will be performed" << endreq;
    }
    else 
    {
        IOverlayDataSvc* overlayInputSvc = SmartIF<IOverlayDataSvc>(tmpService);
        //IOverlayDataSvc* overlayInputSvc = SmartIF<IOverlayDataSvc>(IID_IOverlayDataSvc, tmpService);

        if (!overlayInputSvc) log << MSG::INFO << "Could not find output data service interface" << endreq;
        else                  log << MSG::INFO << "Output data service successfully retrieved" << endreq;
    }

    return sc;
}


StatusCode OverlayTestAlg::execute() 
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), EventModel::EventHeader);

    // For the test job it is necessary to set the time
    eventHeader->setTime(0);

    // Look up the EventOverlay object in the TDS
    SmartDataPtr<Event::EventOverlay> event(eventSvc(), OverlayEventModel::Overlay::EventOverlay);
    if(!event)
    {
        log << MSG::INFO << "Could not find a EventOverlay object in the TDS" << endreq;
    }

    // First, recover CalOverlay objects, to see if we have anything to do
    SmartDataPtr<Event::CalOverlayCol> calOverlayCol(eventSvc(), OverlayEventModel::Overlay::CalOverlayCol);
    if(!calOverlayCol)
    {
        log << MSG::INFO << "Could not find a CalOverlay object in the TDS" << endreq;
    }

    // Now recover any TkrOverlay objects, to see if we have anything to do
    SmartDataPtr<Event::TkrOverlayCol> tkrOverlayCol(eventSvc(), OverlayEventModel::Overlay::TkrOverlayCol);
    if(!tkrOverlayCol)
    {
        log << MSG::INFO << "Could not find a CalOverlay object in the TDS" << endreq;
    }

    // Finally, recover any AcdOverlay objects, to see if we have anything to do
    SmartDataPtr<Event::AcdOverlayCol> acdOverlayCol(eventSvc(), OverlayEventModel::Overlay::AcdOverlayCol);
    if(!acdOverlayCol)
    {
        log << MSG::INFO << "Could not find a CalOverlay object in the TDS" << endreq;
    }

    // Ask service to update event
//    sc = m_overlayInputSvc->selectNextEvent();

    return sc;
}


StatusCode OverlayTestAlg::finalize() 
{    
    return StatusCode::SUCCESS;
}






