// $Header$
#ifndef EVENTCNVSVC_H
#define EVENTCNVSVC_H 1

#include <map>

// Include files
#include "GaudiKernel/ConversionSvc.h"
#include "IEventCnvSvc.h"


// Forward declarations
template <class TYPE> class SvcFactory;


//------------------------------------------------------------------------------
//
// ClassName :    EventCnvSvc
//  
// Description :  Glast event conversion service
//
// Author :       Markus Frank
//
//------------------------------------------------------------------------------


class EventCnvSvc  : virtual public ConversionSvc, virtual public IEventCnvSvc	{

  friend class SvcFactory<EventCnvSvc>;

  /// Map with leaf entries
  typedef std::map<std::string, Leaf*> LeafMap;
  LeafMap m_leaves;


public:

  /// Initialize the service.
  virtual StatusCode initialize();

  /// Finalize the service.
  virtual StatusCode finalize();

  /// Update state of the service
  virtual StatusCode updateServiceState(IOpaqueAddress* pAddress);

  /// Override queryinterface due to enhanced interface
  virtual StatusCode queryInterface(const IID& riid, void** ppvInterface);

  /// Declara data leaf(s)
  virtual StatusCode declareObject(const IEventCnvSvc::Leaf& leaf);

protected:

  /// Standard Constructor
  EventCnvSvc(const std::string& name, ISvcLocator* svc);

  /// Standard Destructor
  virtual ~EventCnvSvc();

};


#endif // EVENTCNVSVC_H
