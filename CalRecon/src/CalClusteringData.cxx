
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "CalClusteringData.h"

CalClusteringData::CalClusteringData( ISvcLocator * svcLocator )
 : m_status(StatusCode::SUCCESS), m_tkrNVertices(0) 
 {
  svcLocator->service("MessageSvc",m_messageSvc,true) ;
  MsgStream log(m_messageSvc,"CalClusteringData::CalClusteringData");

  if (svcLocator->service("EventDataSvc",m_eventSvc,true).isFailure())
   { log<<MSG::ERROR<<"Could not find EventDataSvc"<<endreq ; m_status = StatusCode::FAILURE ; }
  if (svcLocator->service("GlastDetSvc",m_detSvc, true).isFailure())
   { log<<MSG::ERROR<<"Could not find the GlastDetSvc"<<endreq ; m_status = StatusCode::FAILURE ; }

  if (!m_detSvc->getNumericConstByName(std::string("CALnLayer"), &m_calNLayers)) 
   { log<<MSG::ERROR<<"Constant CALnLayer not defined"<<endreq ; m_status = StatusCode::FAILURE ; }
  if (!m_detSvc->getNumericConstByName(std::string("CsIWidth"),&m_calCsIWidth))
   { log<<MSG::ERROR<<"CsIWidth not defined"<<endreq ; m_status = StatusCode::FAILURE ; }
  if (!m_detSvc->getNumericConstByName(std::string("CsIHeight"),&m_calCsIHeight))
   { log<<MSG::ERROR<<"CsIHeight not defined"<<endreq ; m_status = StatusCode::FAILURE ; } 
 }

void CalClusteringData::beginEvent()
 {
  MsgStream log(m_messageSvc,"CalClusteringData::beginEvent") ;

  m_tkrNVertices = 0 ;
  SmartDataPtr<Event::TkrVertexCol> tkrRecData(getEventSvc(),EventModel::TkrRecon::TkrVertexCol) ;
  if (tkrRecData==0)
   { log<<MSG::DEBUG<<"No TKR Reconstruction available "<<endreq ; }
  else
   {
    m_tkrNVertices = tkrRecData->size();
    log<<MSG::DEBUG<<"Number of tracks = "<<m_tkrNVertices<<endreq ;
    if (m_tkrNVertices>0)
     {
      m_tkrFrontVertexDirection = tkrRecData->front()->getDirection() ;
      m_tkrFrontVertexPosition = tkrRecData->front()->getPosition() ;
      m_slope = fabs(m_tkrFrontVertexDirection.z()) ;
      log<<MSG::DEBUG<<"Track direction = "<<m_slope<<endreq ;
     }
    else
     { log<<MSG::DEBUG<<"No reconstructed tracks "<<endreq ; }	
   } 
 }


