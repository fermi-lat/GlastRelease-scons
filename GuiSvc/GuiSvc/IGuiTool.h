#ifndef _H_IGuiTool
#define _H_IGuiTool

// includes
#include "GaudiKernel/IAlgTool.h"

namespace gui {class GuiMgr;}

static const InterfaceID IID_IGuiTool("IGuiTool", 1 , 0); 
 

/** @class IGuiTool
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
 * #include "FluxSvc/IGuiTool.h"
 * #include "FluxSvc/ISpectrumFactory.h"
 * #include "FluxSvc/IFluxSvc.h"
 * // GRB includes
 * 
 * #include "GRB/GRBSpectrum.h"
 * #include "GRBmaker/GRBobsSpectrum.h"
 * 
 * 
 * class RegisterGRB : public AlgTool, virtual public IGuiTool {
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
 *     declareInterface<IGuiTool>(this);
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

class IGuiTool : virtual public IAlgTool {
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IGuiTool; }
    /// Pass a pointer to the gui manager to the tool. 
    virtual StatusCode initialize(gui::GuiMgr*) = 0;
 };

#endif  // _H_IGuiTool
