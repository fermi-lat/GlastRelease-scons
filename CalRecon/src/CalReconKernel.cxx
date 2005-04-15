
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "CalReconKernel.h"

DECLARE_TOOL_FACTORY(CalReconKernel) ;

CalReconKernel::CalReconKernel
 ( const std::string & type, 
   const std::string & name, 
   const IInterface* parent)
 : AlgTool(type,name,parent)
 { declareInterface<CalReconKernel>(this) ; }

CalReconKernel::~CalReconKernel()
 {}

StatusCode CalReconKernel::initialize()
 {
  m_status = StatusCode::SUCCESS ;
  
  //svcLocator->service("MessageSvc",m_messageSvc,true) ;
  service("MessageSvc",m_messageSvc,true) ;
  MsgStream log(m_messageSvc,"CalReconKernel::CalReconKernel");

  if (service("EventDataSvc",m_eventSvc,true).isFailure())
   { log<<MSG::ERROR<<"Could not find EventDataSvc"<<endreq ; m_status = StatusCode::FAILURE ; }
  if (service("GlastDetSvc",m_detSvc, true).isFailure())
   { log<<MSG::ERROR<<"Could not find the GlastDetSvc"<<endreq ; m_status = StatusCode::FAILURE ; }

  if (!m_detSvc->getNumericConstByName(std::string("CALnLayer"), &m_calNLayers)) 
   { log<<MSG::ERROR<<"Constant CALnLayer not defined"<<endreq ; m_status = StatusCode::FAILURE ; }
  else
   {
    // David Chamont: number of layers is also hardcoded in CalCluster !
    //  it will not hurt if I check the two values are consistent.
    Event::CalCluster fake(0,Point()) ;
    if ( fake.getEneLayer().size() != ((unsigned int)m_calNLayers ))
     { log<<MSG::FATAL<<"Inconsistent number of layers"<<endreq ; m_status = StatusCode::FAILURE ; }
   }
  if (!m_detSvc->getNumericConstByName(std::string("CsIWidth"),&m_calCsIWidth))
   { log<<MSG::ERROR<<"CsIWidth not defined"<<endreq ; m_status = StatusCode::FAILURE ; }
  if (!m_detSvc->getNumericConstByName(std::string("CsIHeight"),&m_calCsIHeight))
   { log<<MSG::ERROR<<"CsIHeight not defined"<<endreq ; m_status = StatusCode::FAILURE ; } 

  return m_status ;
 }
 
void CalReconKernel::reviewEvent()
 {
  MsgStream log(m_messageSvc,"CalReconKernel::reviewEvent") ;

  // Ensure CalRecon/Event directory in TDS
  // David Chamont: I tried to move this section in the constructor above,
  //  which is called during the initization of CalClustersAlg, but
  //  for some not investigated reason it does not work.
  DataObject * pnode = 0 ;
  if ((getEventSvc()->retrieveObject(EventModel::CalRecon::Event,pnode)).isFailure()
      && (getEventSvc()->registerObject(EventModel::CalRecon::Event,new DataObject)).isFailure())
   { log<<MSG::ERROR<<"Could not create CalRecon directory in TDS"<<endreq ; m_status = StatusCode::FAILURE ; } 

  // CalXtalRecCol
  m_calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>
   (getEventSvc(),EventModel::CalRecon::CalXtalRecCol) ; 
  if (!m_calXtalRecCol)
   { log<<MSG::VERBOSE<<"No CalXtalRecCol"<<endreq ; }
        
  // CalClusterCol
  m_calClusterCol = SmartDataPtr<Event::CalClusterCol>
   (getEventSvc(),EventModel::CalRecon::CalClusterCol) ;
  if (!m_calClusterCol)
   {
    m_calClusterCol = new Event::CalClusterCol() ;
    if ((getEventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,m_calClusterCol)).isFailure())
     {
      log<<MSG::ERROR<<"Cannot register CalClusterCol"<<endreq ;
      m_status = StatusCode::FAILURE ;
     }
   }
    
  // Tracker
  m_tkrNVertices = 0 ;
  m_tkrSlope = 0 ;
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
      m_tkrSlope = fabs(m_tkrFrontVertexDirection.z()) ;
      log<<MSG::DEBUG<<"Track direction = "<<m_tkrSlope<<endreq ;
     }
    else
     { log<<MSG::DEBUG<<"No reconstructed tracks "<<endreq ; }	
   } 
 }


