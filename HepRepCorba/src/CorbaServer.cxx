#include "CorbaServer.h"
#include "ace/Signal.h"

#include "HepEventServer_impl.h"

#include <fstream>
#include <iostream>

extern int done = 1;

// A signal handler that will close the server object
extern "C" void handler (int)
{
  done = 0;
}

void CorbaServer::run()
{  
  int orb_argc = 4;
  char* orb_argv[] = {"-ORBStdProfileComponents", "0", "-ORBDottedDecimalAddresses", "1"};
  CORBA::ORB_var orb = CORBA::ORB_init(orb_argc, orb_argv);
  // CORBA::ORB_var orb = CORBA::ORB_init(argc,argv);
  // Activate the Portable Object Adapter.
  CORBA::Object_var obj = orb->resolve_initial_references("RootPOA");
  PortableServer::POA_var poa = PortableServer::POA::_narrow(obj);
  PortableServer::POAManager_var mgr = poa->the_POAManager();
  mgr->activate();
    
  // Instantiate the HepEventServer servant.
  HepEventServer_impl* hepEventServer_impl = 
    new HepEventServer_impl(m_svcAdapter, m_registry);  
    
  PortableServer::ObjectId_var oid =
    poa->activate_object( hepEventServer_impl );
    
  CORBA::Object_var hepEventServer_var = 
    poa->id_to_reference( oid.in() );
    
  // Write stringified reference to file hepeventserver.ior.
  CORBA::String_var stringifiedIOR =
    orb->object_to_string(hepEventServer_var.in());

  std::ofstream fout("hepeventserver.ior");
  if (fout.good()) {
    fout << stringifiedIOR.in() << std::endl;
    fout.close( );
  } 

  // hepEventServer_impl->setAppMgr(m_appMgrUI);
  // hepEventServer_impl->setSvcLocator(serviceLocator());
    
  hepEventServer_impl->initHepRep();

  ACE_Sig_Action sa ((ACE_SignalHandler) handler, SIGINT);
    
  while(done)
    {
      // Listen for client requests
      if(orb->work_pending())
	orb->perform_work();
    }
  
  orb->shutdown(1);
  
} 
