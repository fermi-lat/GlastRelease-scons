// $Id$

#ifndef GLASTSVC_HISTORYSVC_H
 #include "GlastSvc/HistorySvc/GlastHistorySvc.h"
#endif

#include "GaudiKernel/HistoryObj.h"

#include "GaudiHistory/JobHistory.h"
#include "GaudiHistory/AlgorithmHistory.h"
#include "GaudiHistory/AlgToolHistory.h"
#include "GaudiHistory/DataHistory.h"
#include "GaudiHistory/ServiceHistory.h"

#ifndef KERNEL_SVCFACTORY_H
 #include "GaudiKernel/SvcFactory.h"
#endif
#ifndef GAUDIKERNEL_ISVCLOCATOR_H
 #include "GaudiKernel/ISvcLocator.h"
#endif
#ifndef GAUDIKERNEL_MSGSTREAM_H
 #include "GaudiKernel/MsgStream.h"
#endif

#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IAlgManager.h"
#include "GaudiKernel/IAlgorithm.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IService.h"

//#include "GaudiKernel/IAlgContextSvc.h"

#include <iostream>
#include <sstream>


static SvcFactory<GlastHistorySvc> s_factory;
const ISvcFactory& GlastHistorySvcFactory = s_factory;

using namespace std;
//extern char** environ;


void env_dummy(char**) {}

GlastHistorySvc::GlastHistorySvc( const std::string& name, ISvcLocator* svc )
  : Service( name, svc ), 
  m_isInitialized(false),
  m_dump(false),
  m_dumpJobHistory(false),
  //p_algCtxSvc(0),
  m_jobHistory(0)
{
    declareProperty("Dump",m_dump);
    declareProperty("DumpJobHistory",m_dumpJobHistory);

    // hack to bring in environ
   // char** tmp = environ;
   // env_dummy(tmp);

}


GlastHistorySvc::~GlastHistorySvc() {
    if (m_jobHistory) {
        delete m_jobHistory;
        m_jobHistory = 0;
    }

    m_algmap.clear();
    m_algtoolmap.clear();
    m_svcmap.clear();
}


StatusCode GlastHistorySvc::initialize() {
  MsgStream log( msgSvc(), name() );

  StatusCode status = Service::initialize();

  log << MSG::DEBUG << "Initializing GlastHistorySvc" << endreq;

    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    status = service ("IncidentSvc", incsvc, true);
    if( status.isFailure() ){
        log << MSG::WARNING << "GlastHistorySvc failed to load IncidentSvc, no property logging" << endreq;
        return StatusCode::SUCCESS;
    }

    incsvc->addListener(this, "BeginEvent", 100);


  // HMK Unused? static const bool CREATEIF(true);

  //if ( service("AlgContextSvc",p_algCtxSvc,CREATEIF).isFailure() ) {
  //  log << MSG::ERROR << "unable to get the AlgContextSvc" << endreq;
  //  return StatusCode::FAILURE;
 // }

  m_jobHistory = new JobHistory;

  /// Get all the Algorithms 


  m_isInitialized = true;

  return status;

}


// handle "incidents"
void GlastHistorySvc::handle(const Incident &inc)
{
    if( inc.type()=="BeginEvent")beginEvent();
}


void GlastHistorySvc::beginEvent() // should be called at the beginning of an event
{ 
    // Cannot necessarily get the collection of algorithms during initialization, as we don't know 
    // the ordering of initialization, so we might miss all or some of the algs, svcs, etc.
    // Also, it seems we cannot wait until finalization to do the dumping of properties, as the algorithms
    // may no longer exist at that point.

    static bool init = false;

    if (init) return;

  MsgStream log( msgSvc(), name() );

  StatusCode sc;
  IAlgManager* algMgr = 0;
  sc = Gaudi::svcLocator()->queryInterface( IID_IAlgManager, 
					    (void**) &algMgr );
  if (sc.isFailure()) {
    log << MSG::WARNING << "Could not get AlgManager" << endreq;
  } else {

    std::list<IAlgorithm*> algs;
    algs = algMgr->getAlgorithms();
    std::list<IAlgorithm*>::const_iterator itr;
    for (itr=algs.begin(); itr!=algs.end(); ++itr) {
      Algorithm* alg = dynamic_cast<Algorithm*> (*itr);
      if (!alg) 
          log << MSG::INFO << "Null Alg" << endreq;
      registerAlg( *alg );
    }
    
    log << MSG::INFO;
    log << "Registered " << algs.size() << " Algorithms" << endreq;
    
  }

  /// Get all the AlgTools

  std::set<const IAlgTool*>::const_iterator itra;
  for (itra = m_ialgtools.begin(); itra != m_ialgtools.end(); ++itra) {
    registerAlgTool(**itra);
  }

  log << MSG::INFO << "Registered " << m_algtools.size() << " AlgTools"
      << endreq;

  /// Get all the Services

  std::list<IService*> svcs = Gaudi::svcLocator()->getServices();
  
  std::list<IService*>::const_iterator itrs;
  for (itrs=svcs.begin(); itrs!=svcs.end(); ++itrs) {
    registerSvc(**itrs);
  }
  
  log << MSG::INFO;
  log << "Registered " << svcs.size() << " Services" << endreq;

 if (m_dump) {
    listProperties();
  }


  init = true;


}



StatusCode GlastHistorySvc::finalize() {
  MsgStream log( msgSvc(), name() );

 // if (m_dump) {
 //   listProperties();
 // }

  if (m_dumpJobHistory) {
    listJobHistory();
  }

  StatusCode status = Service::finalize();

  if ( status.isSuccess() )
    log << MSG::INFO << "Service finalised successfully" << endreq;

  return status;

}


StatusCode GlastHistorySvc::queryInterface(const InterfaceID& riid, 
				      void** ppvInterface) 
{
  StatusCode sc = StatusCode::FAILURE;
  if ( ppvInterface ) {
    *ppvInterface = 0;
    
    if ( IHistorySvc::interfaceID().versionMatch(riid) )    {
      //    if ( riid == IID_HistorySvc ) {
      *ppvInterface = static_cast<IHistorySvc*>(this);
      sc = StatusCode::SUCCESS;
      //      addRef();
    }else if (IID_IIncidentListener.versionMatch(riid) ) {
		*ppvInterface = (IIncidentListener*)this;
        sc = StatusCode::SUCCESS;
    } else
      sc = Service::queryInterface( riid, ppvInterface );    
  }
  return sc;
}


StatusCode GlastHistorySvc::registerAlg(const Algorithm &alg) {

  MsgStream log( msgSvc(), name() );

  m_algs.insert(&alg);
  
  JobHistory *job = getJobHistory();
  AlgorithmHistory *algHist = new AlgorithmHistory(alg, job);
  m_algmap[&alg] = algHist;

  log << MSG::DEBUG << "Registering algorithm: ";
  log.setColor(MSG::CYAN);
  log << alg.name() << endreq;

  return StatusCode::SUCCESS;

}

StatusCode GlastHistorySvc::listProperties(const Algorithm &alg) const {

  MsgStream log( msgSvc(), name() );
  log << MSG::INFO << "Dumping properties for " << alg.name() << endl;

  AlgorithmHistory *hist = getAlgHistory( alg );

  if (hist == 0) {
    return StatusCode::FAILURE;
  }

  ostringstream ost;
  ost << *hist;

  std::string str = ost.str();

  log << MSG::INFO << alg.name() << " --> " << endl << str << endreq;
  
  

  return StatusCode::SUCCESS;

}

AlgorithmHistory* GlastHistorySvc::getAlgHistory(const Algorithm &alg) const {

  MsgStream log( msgSvc(), name() );

  const Algorithm *palg = &alg;
  set<const Algorithm*>::const_iterator itr = m_algs.find(palg);
  if ( itr == m_algs.end() ) {
    log << MSG::ERROR << "Algorithm " << alg.name() << " not registered"
	<< endreq;
    return 0;
  }

  map<const Algorithm*, AlgorithmHistory*>::const_iterator itr2;
  itr2 = m_algmap.find( *itr );

  return ( itr2->second );

}

void
GlastHistorySvc::getAlgHistory(std::set<AlgorithmHistory*>& algs) const {

  set<const Algorithm*>::const_iterator itr;
  for (itr=m_algs.begin(); itr!=m_algs.end(); ++itr) {
    AlgorithmHistory *ah = m_algmap.find(*itr)->second;

    algs.insert(ah);
  }

}


StatusCode GlastHistorySvc::registerJob() {


  return StatusCode::SUCCESS;

}


StatusCode GlastHistorySvc::listProperties() const {

  MsgStream log( msgSvc(), name() );
  log << MSG::INFO;
  log.setColor(MSG::CYAN);
  log << "Dumping properties for all Algorithms" << endreq;

  std::map<const Algorithm*, AlgorithmHistory*>::const_iterator itr;
  for (itr=m_algmap.begin(); itr != m_algmap.end(); ++itr) {
    const Algorithm* alg = itr->first;

    listProperties( *alg );

  }

  log << MSG::INFO;
  log.setColor(MSG::CYAN);
  log << "Dumping properties for all AlgTools" << endreq;

  std::map<const AlgTool*, AlgToolHistory*>::const_iterator itr_a;
  for (itr_a=m_algtoolmap.begin(); itr_a != m_algtoolmap.end(); ++itr_a) {
    const AlgTool* alg = itr_a->first;

    listProperties( *alg );

  }

  log << MSG::INFO;
  log.setColor(MSG::CYAN);
  log << "Dumping properties for all Services" << endreq;

  std::map<const IService*, ServiceHistory*>::const_iterator itr_s;
  for (itr_s=m_svcmap.begin(); itr_s != m_svcmap.end(); ++itr_s) {
    const IService* svc = itr_s->first;

    listProperties( *svc );

  }

  return StatusCode::SUCCESS;
}


StatusCode GlastHistorySvc::listJobHistory() const {

  MsgStream log( msgSvc(), name() );
  log << MSG::INFO;
  log.setColor(MSG::CYAN);
  log << "Dumping properties for Job";
  log.resetColor();
  
  ostringstream ost;
  ost << *m_jobHistory;
  std::string str = ost.str();

  log << std::endl << str << endreq;
  
  return StatusCode::SUCCESS;

}

JobHistory*
GlastHistorySvc::getJobHistory() const {

  return m_jobHistory;

}



DataHistory* GlastHistorySvc::createDataHistoryObj(const CLID& /*id*/, 
                                 const std::string& /*key*/, 
				 const std::string& /*storeName*/) {

                     /*
  MsgStream log( msgSvc(), name() );
  
  log << MSG::DEBUG << "reg DataObj, clid: " << id 
      << "  key: " << key << endreq;

  IAlgorithm* ialg = p_algCtxSvc->currentAlg();
  AlgorithmHistory *algHist;
  if (ialg == 0) {
    log << MSG::DEBUG
	<< "Could not discover current Algorithm:" << endl
	<< "          object CLID: " << id << "  key: \"" << key 
	<< "\"" << endreq;
    algHist = 0;
  } else {

    Algorithm* alg = dynamic_cast<Algorithm*>(ialg);
    if (alg != 0) {
      algHist = getAlgHistory( *alg );
    } else {
      log << MSG::ERROR
	  << "Could not extract concerete Algorithm:"
	  << endl
	  << "          object CLID: " << id << "  key: \"" << key 
	  << "\"" << endreq;
      algHist = 0;
    }
  }

  DataHistory *hist = new DataHistory(id, key, storeName, algHist);

  return hist;
*/
                     return 0;
}



DataHistory* GlastHistorySvc::getDataHistory(const CLID& /*id*/, const std::string& /*key*/, 
			   const std::string& /*storeName*/) const {

  return 0;

}


StatusCode GlastHistorySvc::registerSvc(const IService &svc) {

  MsgStream log( msgSvc(), name() );

  if ( svc.name() == "HistoryStore" ) {
    //    log << MSG::WARNING << "not registering store" << endreq;
    return StatusCode::SUCCESS;
  }

  const IService* psvc = &svc;
  map<const IService*, ServiceHistory*>::const_iterator itr = 
    m_svcmap.find(psvc);
  if (itr == m_svcmap.end()) {

    log << MSG::DEBUG << "Registering Service: ";
    log.setColor(MSG::CYAN);
    log << svc.name() << endreq;

    ServiceHistory *svcHist = new ServiceHistory(&svc);
    m_svcmap[psvc] = svcHist;

  }

  return StatusCode::SUCCESS;

}


ServiceHistory* GlastHistorySvc::getServiceHistory(const IService &svc) const {

  MsgStream log( msgSvc(), name() );

  const IService *psvc = &svc;
  map<const IService*, ServiceHistory*>::const_iterator itr = 
    m_svcmap.find(psvc);
  if ( itr == m_svcmap.end() ) {
    log << MSG::ERROR << "Service " << svc.name() << " not registered"
	<< endreq;
    return 0;
  }

  return ( itr->second );

}

void GlastHistorySvc::getServiceHistory(std::set<ServiceHistory*>& svcs) const {

  set<const IService*>::const_iterator itr;
  for (itr=m_svcs.begin(); itr!=m_svcs.end(); ++itr) {
    ServiceHistory *sh = m_svcmap.find( *itr )->second;

    svcs.insert(sh);
  }

}


StatusCode GlastHistorySvc::listProperties(const IService &svc) const {

  MsgStream log( msgSvc(), name() );
  log << MSG::INFO << "Dumping properties for " << svc.name() << endl;

  ServiceHistory *hist = getServiceHistory( svc );

  if (hist == 0) {
    return StatusCode::FAILURE;
  }

  ostringstream ost;
  ost << *hist;

  std::string str = ost.str();

  log << MSG::INFO << svc.name() << " --> " << endl << str << endreq;
  
  

  return StatusCode::SUCCESS;

}


StatusCode GlastHistorySvc::registerAlgTool(const IAlgTool& ialg) {

  MsgStream log( msgSvc(), name() );

  /*
  if (! m_isInitialized) {
    if (p_algCtxSvc == 0) {
      if ( service("AlgContextSvc",p_algCtxSvc,true).isFailure() ) {
	log << MSG::ERROR << "unable to get the AlgContextSvc" << endreq;
	return StatusCode::FAILURE;
      }
    }

    m_ialgtools.insert(&ialg);
    return StatusCode::SUCCESS;
  }
  */

  const AlgTool *alg = dynamic_cast<const AlgTool*>( &ialg );
  if ( alg == 0 ) {
    log << MSG::ERROR << "Could not dcast IAlgTool \"" << ialg.name()
	<< "\" to an AlgTool" << endreq;
    return StatusCode::FAILURE;
  }

  m_algtools.insert(alg);

  const JobHistory *job = getJobHistory();
  AlgToolHistory *algHist = new AlgToolHistory(*alg, job);
  m_algtoolmap[alg] = algHist;

  log << MSG::DEBUG << "Registering algtool: ";
  log.setColor(MSG::CYAN);
  log << alg->name() << endreq;

  return StatusCode::SUCCESS;

}



StatusCode GlastHistorySvc::listProperties(const IAlgTool& alg) const {

  MsgStream log( msgSvc(), name() );
  log << MSG::INFO << "Dumping properties for " << alg.name() << endl;

  AlgToolHistory *hist = getAlgToolHistory( alg );

  if (hist == 0) {
    return StatusCode::FAILURE;
  }

  ostringstream ost;
  ost << *hist;

  std::string str = ost.str();

  log << MSG::INFO << alg.name() << " --> " << endl << str << endreq;
  
 
  return StatusCode::SUCCESS;

}


AlgToolHistory* GlastHistorySvc::getAlgToolHistory(const IAlgTool& alg) const {
  MsgStream log( msgSvc(), name() );

  const AlgTool *palg = dynamic_cast<const AlgTool*>(&alg);
  set<const AlgTool*>::const_iterator itr = m_algtools.find(palg);
  if ( itr == m_algtools.end() ) {
    log << MSG::ERROR << "AlgTool " << alg.name() << " not registered"
	<< endreq;
    return 0;
  }

  map<const AlgTool*, AlgToolHistory*>::const_iterator itr2;
  itr2 = m_algtoolmap.find( *itr );

  return ( itr2->second );

}


void GlastHistorySvc::getAlgToolHistory(std::set<AlgToolHistory*>& algs) const {

  set<const AlgTool*>::const_iterator itr;
  for (itr=m_algtools.begin(); itr!=m_algtools.end(); ++itr) {
    AlgToolHistory *ah = m_algtoolmap.find(*itr)->second;

    algs.insert(ah);
  }

}
