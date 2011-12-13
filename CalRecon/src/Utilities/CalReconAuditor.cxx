
#include "GaudiKernel/Auditor.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Chrono.h"
#include "GaudiKernel/IChronoStatSvc.h"
#include "GaudiKernel/IAlgorithm.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include <map>

/**   
* @class CalReconAuditor
*
* Experiment the monitoring of algorithms.
*
* $Header$
*/

class CalReconAuditor : virtual public Auditor
 {
  public:

    CalReconAuditor( const std::string &, ISvcLocator * ) ;
    virtual StatusCode initialize() ;
    virtual StatusCode finalize() ;
    virtual ~CalReconAuditor() ;

    //virtual StatusCode beforeInitialize( IAlgorithm * ) ;
    //virtual StatusCode afterInitialize( IAlgorithm * ) ;
    virtual StatusCode beforeExecute( IAlgorithm * ) ;
    virtual StatusCode afterExecute( IAlgorithm * ) ;
    //virtual StatusCode beforeFinalize( IAlgorithm * ) ;
    //virtual StatusCode afterFinalize( IAlgorithm * ) ;

  private :
    
    //! names of algorithms to audit
    StringArrayProperty  m_algoNames ;

    // some pointers to services  
    MsgStream * m_log ;
    IChronoStatSvc * m_chronoSvc ;

    // internals
    std::map< IAlgorithm *, bool > m_algoStatus ;

 } ;


//==========================================================
// construction
//==========================================================

DECLARE_AUDITOR_FACTORY(CalReconAuditor) ;

CalReconAuditor::CalReconAuditor( const std::string & name, ISvcLocator * svcLocator )
  : Auditor(name,svcLocator)
 {
  std::vector<std::string> algoNamesVec ;
  algoNamesVec.push_back("CalClustersAlg") ;
  algoNamesVec.push_back("CalEventEnergyAlg") ;
  declareProperty("algoNames",m_algoNames=algoNamesVec) ;
 }

StatusCode CalReconAuditor::initialize()
 {
  Auditor::initialize() ;
  
  // message stream
  IMessageSvc * messageSvc ;
  service("MessageSvc",messageSvc,true) ;
  m_log = new MsgStream(messageSvc,name()) ;

  // chrono stat svc
  if (service("ChronoStatSvc",m_chronoSvc,true).isFailure())
   {
    (*m_log)<<MSG::ERROR<<"Could not find CalReconSvc"<<endreq ;
    return StatusCode::FAILURE ;
   }

  (*m_log)<<MSG::DEBUG<<"initialize() "<<endreq ;
  return StatusCode::SUCCESS ;
 }


//==========================================================
// events loop
//==========================================================

//StatusCode CalReconAuditor::beforeInitialize( IAlgorithm * alg )
// {
//  (*m_log)<<MSG::DEBUG<<"beforeInitialize("<<alg->name()<<") "<<endreq ;
//  return StatusCode::SUCCESS ;
// }
//
//
//StatusCode CalReconAuditor::afterInitialize( IAlgorithm * algo )
// {
//  (*m_log)<<MSG::DEBUG<<"afterInitialize("<<alg->name()<<") "<<endreq ;
//  return StatusCode::SUCCESS ;
// }
//

StatusCode CalReconAuditor::beforeExecute( IAlgorithm * algo )
 {
  // upgrade m_algoStatus
  std::map< IAlgorithm *, bool >::iterator itr ;
  itr = m_algoStatus.find(algo) ;
  if ( itr == m_algoStatus.end() )
   {
    bool found = false ;
    const std::vector< std::string > & algoNames = m_algoNames ;
    std::vector< std::string >::const_iterator algoName ;
    for ( algoName = algoNames.begin(); algoName != algoNames.end() ; algoName++ ) 
     {
      if ((*algoName)==algo->name())
       {
        found = true ;
        m_algoStatus[algo] = true ;
       }
     }
    if (!found)
     { m_algoStatus[algo] = false ; }
   }

  // look if the algo is under monitoring
  if (m_algoStatus[algo]==false)
   { return StatusCode::SUCCESS ; }

  // real work
  m_chronoSvc->chronoStart(algo->name()) ;
  return StatusCode::SUCCESS ;
 }


StatusCode CalReconAuditor::afterExecute( IAlgorithm * algo )
 {
  // look if the algo is under monitoring
  if (m_algoStatus[algo]==false)
   { return StatusCode::SUCCESS ; }

  // stop chrono
  m_chronoSvc->chronoStop(algo->name()) ;

  // display the last time intervall
  IChronoStatSvc::ChronoTime delta
   = m_chronoSvc->chronoDelta(algo->name(),IChronoStatSvc::USER)/1000 ;
  (*m_log)<<MSG::DEBUG<<algo->name()<<" user time : "<<delta<<" ms"<<endreq ;

  // end
  return StatusCode::SUCCESS ;
 }


//StatusCode CalReconAuditor::beforeFinalize( IAlgorithm * alg )
// {
//  (*m_log)<<MSG::DEBUG<<"beforeFinalize("<<alg->name()<<") "<<endreq ;
//  return StatusCode::SUCCESS ;
// }
//
//
//StatusCode CalReconAuditor::afterFinalize( IAlgorithm * alg )
// {
//  (*m_log)<<MSG::DEBUG<<"afterFinalize("<<alg->name()<<") "<<endreq ;
//  return StatusCode::SUCCESS ;
// }
//

//==========================================================
// destruction
//==========================================================

StatusCode CalReconAuditor::finalize()
 {
  (*m_log)<<MSG::DEBUG<<"finalize() "<<endreq ;
  return Auditor::finalize() ;
 }

CalReconAuditor::~CalReconAuditor()
 {}



