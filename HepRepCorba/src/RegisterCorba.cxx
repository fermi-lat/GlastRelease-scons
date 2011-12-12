/** @file UserSource.cxx 

*/

#include "HepRepSvc/IHepRepSvc.h"
#include "HepRepSvc/IRegister.h"

#include "CorbaServer.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

/** @class Register
*  @brief Registe
*  
*   @author Riccardo Giannitrapani
*   $Header$
*/
class RegisterCorba : public AlgTool, virtual public IRegister{
public:

    RegisterCorba( const std::string& type, 
       const std::string& name, const IInterface* parent);
    virtual ~RegisterCorba() { }

    /// implement to define sources: will be called from FluxSvc
    StatusCode registerMe(IHepRepSvc* );

};


// Static factory for instantiation of algtool objects
//static ToolFactory<RegisterCorba> s_factory;
//const IToolFactory& RegisterCorbaFactory = s_factory;
DECLARE_TOOL_FACTORY(RegisterCorba);

// Standard Constructor
RegisterCorba::RegisterCorba(const std::string& type, 
           const std::string& name, 
           const IInterface* parent)
  : AlgTool( type, name, parent ) 
{
  
  // Declare additional interface
  declareInterface<IRegister>(this);
}


StatusCode RegisterCorba::registerMe(IHepRepSvc* hrsvc) 
{
  MsgStream  log(msgSvc(), name());
  log << MSG::INFO << "Register corba server ..." << endreq;
  
  CorbaServer* corba = new CorbaServer();
  hrsvc->setServer(corba);
  
  return StatusCode::SUCCESS;
} 
