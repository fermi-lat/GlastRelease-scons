// $Id$
// Stolen from GaudiSvc in Gaudi v18r1, as it was inactive by default and we may desire
// to customize this version to allow the properties to be stored in ROOT

#ifndef GLASTSVC_HISTORYSVC_H
#define GLASTSVC_HISTORYSVC_H

#include "GaudiKernel/IHistorySvc.h"

#ifndef GAUDIKERNEL_SERVICE_H
 #include "GaudiKernel/Service.h"
#endif
#ifndef KERNEL_STATUSCODES_H
 #include "GaudiKernel/StatusCode.h"
#endif

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgTool.h"

#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"


#include <vector>
#include <string>
#include <set>
#include <map>
#include "GaudiKernel/ClassID.h"

class StoreGateSvc;
//class IAlgContextSvc;
class HistoryObj;
//class DataHistory;
class AlgorithmHistory;
class ServiceHistory;

// Forward declarations
template <class TYPE> class SvcFactory;

/** @class GlastHistorySvc GlastHistorySvc.h
 * 
 *  GlastHistorySvc class definition
 *
 */

class GlastHistorySvc: public Service, 
		  virtual public IHistorySvc,
          virtual public IIncidentListener {

public:

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  // Query the interfaces.
  virtual StatusCode queryInterface( const InterfaceID& riid, 
				     void** ppvInterface );
  /// Handles incidents, implementing IIncidentListener interface
  virtual void handle(const Incident& inc);    
	


  virtual StatusCode registerJob();
  virtual StatusCode listProperties() const;
  virtual StatusCode listJobHistory() const;
  virtual JobHistory* getJobHistory() const;

  virtual StatusCode registerSvc(const IService&);
  virtual StatusCode listProperties(const IService&) const;
  virtual ServiceHistory* getServiceHistory(const IService&) const;
  virtual void getServiceHistory(std::set<ServiceHistory*>&) const;

  virtual StatusCode registerAlg(const Algorithm&);
  virtual StatusCode listProperties(const Algorithm&) const;
  virtual AlgorithmHistory* getAlgHistory(const Algorithm&) const;
  virtual void getAlgHistory(std::set<AlgorithmHistory*>&) const;

  virtual StatusCode registerAlgTool(const IAlgTool&);
  virtual StatusCode listProperties(const IAlgTool&) const;
  virtual AlgToolHistory* getAlgToolHistory(const IAlgTool&) const;
  virtual void getAlgToolHistory(std::set<AlgToolHistory*>&) const;

  
  virtual DataHistory* createDataHistoryObj(const CLID& id, 
					    const std::string& key,
					    const std::string& store);
  virtual DataHistory* getDataHistory(const CLID& id, const std::string& key,
				      const std::string& store) const;



protected:

  GlastHistorySvc( const std::string& name, ISvcLocator* svc );

  // Destructor.
  virtual ~GlastHistorySvc();

private:

  void beginEvent();

  // Allow SvcFactory to instantiate the service.
  friend class SvcFactory<GlastHistorySvc>;

  bool m_isInitialized;
  bool m_dump;
  bool m_dumpJobHistory;

  //IAlgContextSvc *p_algCtxSvc;

  std::set<const Algorithm*> m_algs;
  std::map<const Algorithm*, AlgorithmHistory*> m_algmap;

  std::set<const IAlgTool*> m_ialgtools;
  std::set<const AlgTool*> m_algtools;
  std::map<const AlgTool*, AlgToolHistory*> m_algtoolmap;

  std::set<const IService*> m_svcs;
  std::map<const IService*, ServiceHistory*> m_svcmap;

  JobHistory *m_jobHistory;

};

#endif
