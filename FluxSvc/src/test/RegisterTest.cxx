
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/IFluxSvc.h"

class RegisterTest : public AlgTool, virtual public IRegisterSource {
public:
    
    RegisterTest( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RegisterTest() { }
    
    /// implement to define sources: will be called from FluxSvc
    StatusCode registerMe(IFluxSvc* );
    
};


// Static factory for instantiation of algtool objects
static ToolFactory<RegisterTest> s_factory;
const IToolFactory& RegisterTestFactory = s_factory;

// Standard Constructor
RegisterTest::RegisterTest(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool( type, name, parent ) {
    
    // Declare additional interface
    declareInterface<IRegisterSource>(this);
        
}


StatusCode RegisterTest::registerMe(IFluxSvc* fsvc) 
{
    MsgStream  log(msgSvc(), name());
    log << MSG::INFO << "Register test called..." << endreq;
/* the following is an example
    static RemoteSpectrumFactory<GRBSpectrum> factory(fsvc);
*/
    log << MSG::INFO << "About to return" << endreq;
    return StatusCode::SUCCESS;
} 
