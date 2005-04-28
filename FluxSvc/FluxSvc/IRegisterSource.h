#ifndef _H_IRegisterSource
#define _H_IRegisterSource

// includes
#include "GaudiKernel/IAlgTool.h"
#include "flux/ISpectrumFactory.h"
#include <typeinfo>


class IFluxSvc;

static const InterfaceID IID_IRegisterSource("IRegisterSource", 1 , 0); 


/** @class IRegisterSource
* @brief Abstract definition of a tool to be called from FluxSvc to load external ISpectrumFactory enteries
* 
* 
* $Header$
* 
* <br> Example of an implementation:
*   <pre>
*   
* #include "GaudiKernel/AlgTool.h"
* #include "GaudiKernel/MsgStream.h"
* #include "GaudiKernel/ToolFactory.h"
* 
* #include "FluxSvc/IRegisterSource.h"
* #include "FluxSvc/ISpectrumFactory.h"
* #include "FluxSvc/IFluxSvc.h"
* // GRB includes
* 
* #include "GRB/GRBSpectrum.h"
* #include "GRBmaker/GRBobsSpectrum.h"
* 
* 
* class RegisterGRB : public AlgTool, virtual public IRegisterSource {
*  public:
*      
*    RegisterGRB( const std::string& type, const std::string& name, const IInterface* parent);
*     virtual ~RegisterGRB() { }
*     
*     /// implement to define sources: will be called from FluxSvc
*     StatusCode registerMe(IFluxSvc* );
*     
* };
* 
* 
* // Static factory for instantiation of algtool objects
* static ToolFactory<RegisterGRB> s_factory;
* const IToolFactory& RegisterGRBFactory = s_factory;
* 
* // Standard Constructor
* RegisterGRB::RegisterGRB(const std::string& type, 
*                          const std::string& name, 
*                          const IInterface* parent)
*                          : AlgTool( type, name, parent ) {
*     
*     // Declare additional interface
*     declareInterface<IRegisterSource>(this);
*         
* }
* 
* 
* StatusCode RegisterGRB::registerMe(IFluxSvc* fsvc) 
* {
*     MsgStream  log(msgSvc(), name());
*     log << MSG::INFO << "Register GRB Spectra..." << endreq;
*     static RemoteSpectrumFactory<GRBSpectrum> factory(fsvc);
*     
*     log << MSG::INFO << "Register  Sandhia GRB Spectra..." << endreq;
*     static RemoteSpectrumFactory<GRBobsSpectrum> factory2(fsvc);
*     
*     return StatusCode::SUCCESS;
* } 
*</pre>
*/

class IRegisterSource : virtual public IAlgTool {
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IRegisterSource; }
    /// Pass a pointer to the service to the tool. 
    virtual StatusCode registerMe(IFluxSvc*) = 0;
};

template <class T> class RemoteSpectrumFactory : public ISpectrumFactory 
{
public:
    
    RemoteSpectrumFactory(IFluxSvc* svc){
        //Get class name using RTTI:
        std::string m_classname = typeid(T).name();
        int s = m_classname.find_first_of("class");
        if( s==0 ) s=6; //found it
        else s =m_classname.find_first_not_of("0123456789");
        m_classname = m_classname.substr(s);
        svc->addFactory(m_classname, this); 
    }
    /*! return a new Spectrum object
    @param params String from the xml
    @param engine random engine to use
    */

    virtual ISpectrum* instantiate(const std::string& params) const{
        return new T(params);
    }
    
    //! dummy to follow Gaudi model
    virtual void addRef()const{}

    std::string name()const{return m_classname;}
private:
    std::string m_classname;
};


#endif  // _H_IRegisterSource
