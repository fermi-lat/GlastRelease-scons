
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/IFluxSvc.h"
#include "CLHEP/Random/Random.h"

/*! @class RegisterTest
    @brief Example of special tool used to register sources
*/

class RegisterTest : public AlgTool, virtual public IRegisterSource {
public:
    
    RegisterTest( const std::string& type, const std::string& name, const IInterface* parent);
   
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
/*! @class TestSpectrum

*/

class TestSpectrum  : public ISpectrum 
{
    
public:

    /*! This initializes the simulation parsing the parameters.
    
    \param params are set in the xml source library in xml directory.
    An example the xml source declaration for this spectrum should appears:
    \verbatim
    <source name="TestSpectrum">
       <spectrum escale="GeV"> 
          <SpectrumClass name="testSpectrum" params="50 100"/>
          <use_spectrum/>
       </spectrum> 
    </source>
    \endverbatim
  */
   TestSpectrum(const std::string& params)
   {
       std::cout << "TestSpectrum called with parameter string " << params << std::endl;
   }

    ///  particle name that must be known to the particle service
    virtual const char * particleName()const {return "mu+";}
    
    /// calculate the flux, particles/m^2/sr. (default zero)
    virtual double    flux (double time ) const{return 0;}
    
    /// return effective solid angle that will be used to determine the actual rate 
    virtual double solidAngle()const{return  1.0;} //flag that doesen't calculate.
    
    /// return a title describing the spectrum	
    virtual std::string title()const{ return std::string("TestSpectrum"); }
    
    /// a (randomized) interval to the next event.  default is exponential( 1/rate() )
    /// needs to know the cross-sectional area?
    virtual double interval (double time){return 1.;}
    
    /// interface for energy and direction (
    virtual double energy( double time=0){return 2.;}

    virtual std::pair<double,double> dir(double energy){
        return std::make_pair(0,0);
    }
    
};
StatusCode RegisterTest::registerMe(IFluxSvc* fsvc) 
{
    MsgStream  log(msgSvc(), name());
    log << MSG::INFO << "Register test called..." << endreq;
    HepRandomEngine* fluxengine = fsvc->getRandomEngine();
    log << MSG::DEBUG << "My random engine is " << HepRandom::getTheEngine() << endreq;
    log << MSG::DEBUG << "The FluxSvc random engine is " << fluxengine << endreq;

    // this is what should be done, so that the local pointer to the engine is the same as FluxSvc
    HepRandom::setTheEngine(fluxengine);

    static RemoteSpectrumFactory<TestSpectrum> factory(fsvc);
    return StatusCode::SUCCESS;
} 
