#ifndef HepEventServer_impl_hh
#define HepEventServer_impl_hh

#include "HepEventServer.hh"
#include "HepRep.hh"
#include "HepRepSvc/ISvcAdapter.h"
#include "HepRepSvc/IRegistry.h"

class HepRep_impl;

class HepEventServer_impl : public POA_HepEventServer, 
                            public PortableServer::RefCountServantBase
{

public:
  HepEventServer_impl (ISvcAdapter*, IRegistry*);
  ~HepEventServer_impl();
  
  virtual char* attach(const char* clientDesc);
  
  void release();
  
  char* setEvent (const char* eventID);

  CORBA::Long getNumberOfEvents();

  char* getEventTitle();

  /*  
  CORBA::Any* getEventData(const char* dataFormat);
  */
  
  HepRep_ptr getEventData(const char* dataFormat);

  void initHepRep();

  HepRep_impl* hepRep_impl;

protected:
  ISvcAdapter* m_svcAdapter;
  IRegistry* m_registry;
  std::string m_eventID;
  HepRep_var _hepRep_var;
  bool _busy;
  bool _singleUse;
  char* _defaultEventID;
  char* _clientDesc;
};
#endif
