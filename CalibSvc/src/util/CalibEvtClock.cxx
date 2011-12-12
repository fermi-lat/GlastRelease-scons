//$Header$
#include <stdio.h>

#include "CalibEvtClock.h"
#include "CalibData/CalibTime.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"


//  #include "GaudiKernel/TimePoint.h"

/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<CalibEvtClock> Factory;
//const IAlgFactory& CalibEvtClockFactory = Factory;
DECLARE_ALGORITHM_FACTORY(CalibEvtClock);

CalibEvtClock::CalibEvtClock( const std::string&  name, 
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


StatusCode CalibEvtClock::initialize() {
  StatusCode sc;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initialize()" << endreq;

  IDataProviderSvc* calibSvc;
  sc = service("CalibDataSvc", calibSvc, true);

  // Query the IDetDataSvc interface of the calib data service
  sc = calibSvc->queryInterface(IDetDataSvc::interfaceID(), 
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


StatusCode CalibEvtClock::execute( ) {

  MsgStream log(msgSvc(), name());

  // Increment the event counter
  m_eventNumber++;
  log << MSG::INFO << "CalibEvtClock Event number: " << m_eventNumber << endreq;

  // Set the event time
  facilities::Timestamp time = i_evtTime();
  log << MSG::INFO << "Event time: "
      << time.getString()
      << endreq; 
    //      << " Julian day number "
    //      << time.getJulian()
  CalibData::CalibTime ctime(time);
  log << MSG::INFO << "Event time (hours) " << ctime.hours() << endreq;
  //  m_detDataSvc->setEventTime(CalibData::CalibTime(time));
  m_detDataSvc->setEventTime(ctime.getGaudiTime());

  return StatusCode::SUCCESS;
}



StatusCode CalibEvtClock::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "------------- FINALIZE!! -------------------------------------------"
      << endreq;
  log << MSG::INFO << "Total #events: " << m_eventNumber << endreq;
  
  return StatusCode::SUCCESS;
}

facilities::Timestamp CalibEvtClock::i_evtTime( ) {

  return facilities::Timestamp(m_startTime + 
                               (m_eventNumber - 1) * m_delayTime);

}
