/// An HepEvent server; mainly derived from Joseph Perl one .. adapted for GLAST and
/// GAUDI framework

#include "HepRep_impl.hh"
#include "HepEventServer_impl.h"
#include <iostream>
#include <string>

#include <sstream>

HepEventServer_impl::HepEventServer_impl (ISvcAdapter* sa, IRegistry* r)
{
  m_svcAdapter = sa;
  m_registry = r;
  _busy = false;
  // Parse arguments.
  _singleUse = false;
  _defaultEventID = new char[1];
  strcpy(_defaultEventID,"");
}


HepEventServer_impl::~HepEventServer_impl ()
{

}


void HepEventServer_impl::initHepRep()
{

  // Instantiate the HepRep CORBA implementation.
  hepRep_impl = new HepRep_impl();
  _hepRep_var = hepRep_impl -> _this();

  hepRep_impl->setRegistry(m_registry);

}


char* HepEventServer_impl::attach(const char* clientDesc)
{
  if (_busy) {
    return CORBA::string_dup("server unavailable: busy with another client");
  } else {
    _busy = true;
    
    _clientDesc = new char[strlen(clientDesc) + 1];
    strcpy(_clientDesc, clientDesc);

    int descLength = 96;

    //    char* serverDesc = new char[descLength];
    std::string serverDesc;
    serverDesc = 
      "GlastServer version=0.1\ndataFormat=HepRep\nsetEventCommand=next,fluxes,source\n";

    
    std::cout << "Wow, I've been requested by " 
              << _clientDesc << std::endl;

    return CORBA::string_dup(serverDesc.c_str());
  }
}


void HepEventServer_impl::release()
{
  _busy = false;
  std::cout << "I've been abandoned" << std::endl;
}


/// Set the event to the command 
char* HepEventServer_impl::setEvent(const char* command)
{
  // This method is used by now to send command to the server (since in the
  // original implementation of the HepEventServer interface there is not such a
  // method, we decided with Joseph to use setEvent). It accept a string and
  // return a string. For now 4 commands are supported: 1) next -> step to the
  // next event 2) fluxes -> return a list of known fluxes (the list names are
  // separated by ";") 3) source:xxx -> set the FluxSource to xxx 4) eventName
  // -> return the actual event name (for now this name is just Event-x, with x
  // an increasing number; this will be changed in future release).
  std::stringstream sName;
  static temp = 0;
  std::string cmd(command);
  
  std::cout << "I've received from " << _clientDesc << ", setEvent to " 
            << cmd << std::endl;
  
  static unsigned int i = 0;
  
  std::string nextEventMsg; 
  
  if (cmd == "next")
    {
      // Usare l'heprepsvc
      m_svcAdapter->nextEvent(1);
      sName << "Event-" << temp;
      temp++;

      m_eventID = sName.str();
      nextEventMsg = "Event set to next";      
    }
  /*
  else if (cmd == "fluxes")
    {
#ifdef  DEFECT_NO_STRINGSTREAM
  std::strstream sNames;
#else
  std::stringstream sNames;
#endif
      
      std::list<std::string>::const_iterator s;
      std::list<std::string> names = m_fluxSvc->fluxNames();

      for(s=names.begin(); s!=names.end(); s++)
        sNames << (*s) << ";";
      
	  sNames << '\0';
      nextEventMsg = sNames.str();
    }
  else if (cmd.substr(0,7) == "source:")
    {
      std::string source = cmd.substr(7,cmd.size());
      std::cout << "Set the source to " << source << std::endl;

      IFlux *   iflux;
      m_fluxSvc->source(source, iflux);
      nextEventMsg = "Surce changed";      
    }
  */
  else if (cmd == "eventName")
    {
      nextEventMsg = m_eventID.c_str();
    }
  else
    {
      nextEventMsg = "Nope, don't know the message ..\n";      
    }

  return CORBA::string_dup(nextEventMsg.c_str());
}


/// Return the number of events (if known); disabled for now
int HepEventServer_impl::getNumberOfEvents()
{

  return 0;
}


/// Return the title of the current event
char* HepEventServer_impl::getEventTitle()
{
    return CORBA::string_dup(m_eventID.c_str());
}

/// Return the HepRep object
HepRep_ptr HepEventServer_impl::getEventData (const char* dataFormat)
{
  std::cout << "I've just given to " << _clientDesc << " an HepRep object" << std::endl;

  return HepRep::_duplicate (_hepRep_var);
}
