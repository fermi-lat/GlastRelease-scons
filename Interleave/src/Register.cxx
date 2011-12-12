/** @file Register.cxx

*/

#include "SampledBackground.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/IFluxSvc.h"


#include "CLHEP/Random/Random.h"


/** @class RegisterSampledBackground
 *  @brief Register the SampledBackground
 *  
 *   @author Toby Burnett
 *   $Header$
 */
class RegisterSampledBackground : public AlgTool, virtual public IRegisterSource {
public:

    RegisterSampledBackground( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RegisterSampledBackground() { }


    /// implement to define sources: will be called from FluxSvc
    StatusCode registerMe(IFluxSvc* );
    StatusCode initialize(); // overload

};

// Static factory for instantiation of algtool objects
//static ToolFactory<RegisterSampledBackground> s_factory;
//const IToolFactory& RegisterSampledBackgroundFactory = s_factory;
DECLARE_TOOL_FACTORY(RegisterSampledBackground);

// Standard Constructor
RegisterSampledBackground::RegisterSampledBackground(const std::string& type, 
                                                   const std::string& name, 
                                                   const IInterface* parent)
                                                   : AlgTool( type, name, parent ) 
{
    // Declare additional interface
    declareInterface<IRegisterSource>(this);
}


StatusCode RegisterSampledBackground::registerMe(IFluxSvc* fsvc) 
{
    MsgStream  log(msgSvc(), name());
    log << MSG::INFO << "Register SampledBackground" << endreq;

    //declare the factories here:
    static RemoteSpectrumFactory<SampledBackground> factory(fsvc);

    //  needs to use the same random engine as FluxSvc ?
    CLHEP::HepRandomEngine* engine = fsvc->getRandomEngine();
    CLHEP::HepRandom::setTheEngine(engine);

    // Get the initial location from FluxSvc and store in the CrLocation singleton
    //?CrLocation::instance()->setFluxSvc(fsvc);
    initialize(); //?

    return StatusCode::SUCCESS;
} 
StatusCode RegisterSampledBackground::initialize()
{
    //Set the properties
    setProperties();

    return StatusCode::SUCCESS;
}





