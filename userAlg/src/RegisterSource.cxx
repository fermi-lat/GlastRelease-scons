/** @file UserSource.cxx 

*/
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/IFluxSvc.h"

/** @class UserSource
    @brief a test particle source: must inherit from ISpectrum
    */
class UserSource :  public ISpectrum {
public:
    UserSource(const std::string& params, HepRandomEngine* engine=0){};
    
   
    ///  specify correct particle type
    virtual const char * particleName()const{return "mu+";}
    
    
    /// return a title describing the spectrum	
    virtual std::string title()const{return "test spectrum";}
    
    /// a randomized interval to the next event - default is 1/rate()
    virtual double interval (double time) { return 1.0;}
    
    
    /// use the engine and the current time to return an energy
    virtual double energySrc(HepRandomEngine* engine, double time){return 2.;};// in GeV for now

    /// given the selected energy, return a pair(costheta, phi)
    virtual std::pair<double,double> dir(double energy, HepRandomEngine* engine){
        return dir(energy);
    }
     /// calculate the flux, particles/m^2/sr. (default zero)
    virtual double    flux (double time ) const{
        return 0;
    }

    virtual std::pair<float,float> dir(float energy)const {
        return std::make_pair(-1.0,0) ;}

    // not really needed by the interface
    virtual float operator()(float /*r*/)const{return 0;}
    

private:
};


/** @class RegisterSource
 *  @brief Register a sources
 *  
 *   @author Toby Burnett
 *   $Header$
 */
class RegisterSource : public AlgTool, virtual public IRegisterSource {
public:
    
    RegisterSource( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RegisterSource() { }
    
    /// implement to define sources: will be called from FluxSvc
    StatusCode registerMe(IFluxSvc* );
    
};


// Static factory for instantiation of algtool objects
static ToolFactory<RegisterSource> s_factory;
const IToolFactory& RegisterSourceFactory = s_factory;

// Standard Constructor
RegisterSource::RegisterSource(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool( type, name, parent ) {
    
    // Declare additional interface
    declareInterface<IRegisterSource>(this);
        
}


StatusCode RegisterSource::registerMe(IFluxSvc* fsvc) 
{
    MsgStream  log(msgSvc(), name());
    log << MSG::INFO << "Register test Spectra..." << endreq;
    static RemoteSpectrumFactory<UserSource> factory(fsvc);
    
    return StatusCode::SUCCESS;
} 
