//$Header$
#include <stdio.h>

#include "EvtClock.h"
#include "CalibData/CalibTime.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"


//  #include "GaudiKernel/TimePoint.h"

/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<EvtClock> Factory;
const IAlgFactory& EvtClockFactory = Factory;


EvtClock::EvtClock( const std::string&  name, 
		    ISvcLocator*        pSvcLocator )
  : Algorithm     ( name, pSvcLocator )
  , m_eventNumber ( 0 )
  , m_detDataSvc  ( 0 )
{
  declareProperty( "startTime",  
                   m_startTimeAsc = "2003-1-10_00:20");

  // = facilities::Timestamp("2003-1-11").getClibTime());
  declareProperty( "delayTime",  m_delayTime = 2000);
}


StatusCode EvtClock::initialize() {
  StatusCode sc;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initialize()" << endreq;

  IDataProviderSvc* calibSvc;
  sc = service("CalibDataSvc", calibSvc, true);

  // Query the IDetDataSvc interface of the calib data service
  sc = calibSvc->queryInterface(IID_IDetDataSvc, 
                                (void**) &m_detDataSvc);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not query IDetDataSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  } else {
    log << MSG::DEBUG 
	<< "Retrieved IDetDataSvc interface of CalibDataSvc" 
	<< endreq;
  }

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  unsigned int underpos = m_startTimeAsc.find("_");
  if (underpos < m_startTimeAsc.size()) {
    m_startTimeAsc.replace(underpos, 1, " ");
  }
  m_startTime = facilities::Timestamp(m_startTimeAsc).getClibTime();

  log << MSG::DEBUG << "Properties were read from jobOptions" << endreq;
  log << MSG::INFO << "Time of first event: (ascii) "
      << m_startTimeAsc       << endreq; 
  log << MSG::INFO << "Time of first event: (seconds since 1970) "
      << m_startTime       << endreq; 
  log << MSG::INFO << "Time between events (seconds): "
      << m_delayTime 
      << endreq;
  return StatusCode::SUCCESS;

}


StatusCode EvtClock::execute( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "------------- NEW EVENT! -------------------------------------------"
      << endreq;
  log << MSG::INFO << "Execute()" << endreq;

  // Increment the event counter
  m_eventNumber++;
  log << MSG::INFO << "EvtClock Event number: " << m_eventNumber << endreq;

  // Set the event time
  facilities::Timestamp time = i_evtTime();
  log << MSG::INFO << "Event time: "
      << time.getString()
      << " Julian day number "
      << time.getJulian()
      << endreq; 
  m_detDataSvc->setEventTime(CalibData::CalibTime(time));

  return StatusCode::SUCCESS;
}



StatusCode EvtClock::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "------------- FINALIZE!! -------------------------------------------"
      << endreq;
  log << MSG::INFO << "Total #events: " << m_eventNumber << endreq;
  
  return StatusCode::SUCCESS;
}


//longlong EvtClock::i_evtTime( ) {
facilities::Timestamp EvtClock::i_evtTime( ) {

  return facilities::Timestamp(m_startTime + 
                               (m_eventNumber - 1) * m_delayTime);

}
