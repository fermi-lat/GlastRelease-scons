/** @file UserSource.cxx 

*/

#include "HepRepSvc/IHepRepSvc.h"
#include "HepRepSvc/IRegister.h"

#include "XmlStreamer.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

/** @class Register
*  @brief Registe
*  
*   @author Riccardo Giannitrapani
*   $Header$
*/
class RegisterXml : public AlgTool, virtual public IRegister{
public:

    RegisterXml( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RegisterXml() { }

    /// implement to define sources: will be called from FluxSvc
    StatusCode registerMe(IHepRepSvc* );

};


// Static factory for instantiation of algtool objects
//static ToolFactory<RegisterXml> s_factory;
//const IToolFactory& RegisterXmlFactory = s_factory;
DECLARE_TOOL_FACTORY(RegisterXml);

// Standard Constructor
RegisterXml::RegisterXml(const std::string& type, 
		   const std::string& name, 
		   const IInterface* parent)
                   : AlgTool( type, name, parent ) 
{

    // Declare additional interface
    declareInterface<IRegister>(this);

}


StatusCode RegisterXml::registerMe(IHepRepSvc* hrsvc) 
{
    MsgStream  log(msgSvc(), name());
    log << MSG::INFO << "Register xml streamer ..." << endreq;
    
    XmlStreamer* xml = new XmlStreamer();
    hrsvc->addStreamer("xml",xml);
    
    return StatusCode::SUCCESS;
} 
