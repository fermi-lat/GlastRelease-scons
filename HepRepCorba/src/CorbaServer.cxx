#include "CorbaServer.h"
#include "HepEventServer_impl.h"

#include <fstream>
#include <iostream>

extern int done = 1;


void CorbaServer::run()
{  
  int orb_argc = 1 ;
  char* orb_argv[] = {""};

  const char* options[][2] = { { "giopMaxMsgSize", "40097152" }, { 0, 0 } }; 

  try {
  m_orb = CORBA::ORB_init(orb_argc, orb_argv, "omniORB4",options);
  // Activate the Portable Object Adapter.
  CORBA::Object_var obj = m_orb->resolve_initial_references("RootPOA");
  PortableServer::POA_var poa = PortableServer::POA::_narrow(obj);
    
  // Instantiate the HepEventServer servant.
  HepEventServer_impl* hepEventServer_impl = 
    new HepEventServer_impl(m_svcAdapter, m_registry);  
    
  PortableServer::ObjectId_var oid =
    poa->activate_object( hepEventServer_impl );

  obj = hepEventServer_impl->_this();
  CORBA::String_var sior(m_orb->object_to_string(obj));  

  // I need this to build system command 
  #ifdef WIN32
  std::string separator("\\");
  #else // gcc < 2.9.5.3 dosen't know <sstream>,stringstream
  std::string separator("/");
  #endif
  
  // We build the name of the ior file; if FREDHOME env variable is defined, it
  // means we want to launch FRED automatically after the CORBA initialization;
  // the ior file is than saved in the TMP directory if defined; otherwise the
  // ior file is saved in the actual directory of execution
  std::string iorFileName;  
  if ((::getenv("FREDHOME") != NULL) && (::getenv("TMP") != NULL))
  {
    iorFileName = std::string(::getenv("TEMP")) + separator; 
  }
  iorFileName = iorFileName + "hepeventserver.ior";
      
  std::ofstream fout(iorFileName.c_str());
  if (fout.good()) {
    fout << (char*)sior << std::endl;
    fout.close( );
  }
  
  hepEventServer_impl->_remove_ref();

  PortableServer::POAManager_var mgr = poa->the_POAManager();
  mgr->activate();

  
  hepEventServer_impl->initHepRep();

  // If FREDHOME and TMP are defined, we will try to start FRED automatically
  // passing to it the ior file name
  if ((::getenv("FREDHOME") != NULL) && (::getenv("TMP") != NULL))
  {    
    std::string root("ruby ");
    root = root + std::string(::getenv("FREDHOME"));
    root = root + separator + std::string("fred.rb -s ");
    system((root + iorFileName).c_str());
  }

  // We start the CORBA server
  m_orb->run();    
  m_orb->destroy();
  }
 catch(CORBA::COMM_FAILURE&) {
   std::cerr << "Caught system exception COMM_FAILURE -- unable to contact the "
         << "object." << std::endl;
  }
  catch(CORBA::SystemException&) {
    std::cerr << "Caught CORBA::SystemException." << std::endl;
  }
  catch(CORBA::Exception&) {
    std::cerr << "Caught CORBA::Exception." << std::endl;
  }
  catch(omniORB::fatalException&) {
    std::cerr << "Caught omniORB::fatalException:" << std::endl;
  }
}

void CorbaServer::shutDown()
{
  std::cout << "Destro" << std::endl;
  m_orb->shutdown(0);
}
