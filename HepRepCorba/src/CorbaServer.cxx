#include "CorbaServer.h"

#include "HepEventServer_impl.h"

#include <fstream>
#include <iostream>

extern int done = 1;


void CorbaServer::run()
{  
  int orb_argc = 1 ;
  char* orb_argv[] = {""};

  try {
  CORBA::ORB_var orb = CORBA::ORB_init(orb_argc, orb_argv);
  // Activate the Portable Object Adapter.
  CORBA::Object_var obj = orb->resolve_initial_references("RootPOA");
  PortableServer::POA_var poa = PortableServer::POA::_narrow(obj);
    
  // Instantiate the HepEventServer servant.
  HepEventServer_impl* hepEventServer_impl = 
    new HepEventServer_impl(m_svcAdapter, m_registry);  
    
  PortableServer::ObjectId_var oid =
    poa->activate_object( hepEventServer_impl );

  obj = hepEventServer_impl->_this();
  CORBA::String_var sior(orb->object_to_string(obj));  


  std::ofstream fout("hepeventserver.ior");
  if (fout.good()) {
    fout << (char*)sior << std::endl;
    fout.close( );
  }
  
  hepEventServer_impl->_remove_ref();

  // hepEventServer_impl->setAppMgr(m_appMgrUI);
  // hepEventServer_impl->setSvcLocator(serviceLocator());

  PortableServer::POAManager_var mgr = poa->the_POAManager();
  mgr->activate();

  
  hepEventServer_impl->initHepRep();

  orb->run();    
  orb->destroy();
 }
 catch(CORBA::COMM_FAILURE& ex) {
   std::cerr << "Caught system exception COMM_FAILURE -- unable to contact the "
         << "object." << std::endl;
  }
  catch(CORBA::SystemException&) {
    std::cerr << "Caught CORBA::SystemException." << std::endl;
  }
  catch(CORBA::Exception&) {
    std::cerr << "Caught CORBA::Exception." << std::endl;
  }
  catch(omniORB::fatalException& fe) {
    std::cerr << "Caught omniORB::fatalException:" << std::endl;
  }
} 
