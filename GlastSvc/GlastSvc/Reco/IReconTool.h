
#ifndef _H_IReconTool
#define _H_IReconTool

#include "GaudiKernel/IAlgTool.h"



// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IReconTool("IReconTool", 1 , 0); 

/** @class IReconTool
* @brief Abstract interface for GLAST reconstruction returning a quantity of interest for cuts.
*
* @author Toby Burnett
* $Header$
* <br> Example of an implementation:
*   <pre>
*   
* #include "GaudiKernel/AlgTool.h"
* #include "GaudiKernel/MsgStream.h"
* #include "GaudiKernel/ToolFactory.h"
* 
* #include "GlastSvc/IReconTool.h"
* 
* 
* class MyRecon : public AlgTool, virtual public IReconTool {
*  public:
*      
*    MyRecon( const std::string& type, const std::string& name, const IInterface* parent);
*     virtual ~RegisterGRB() { }
*     
*     /// implement to get something by name
*     virtual StatusCode get(const std::string& name, double& x)
*
*     
* };
* 
* 
* // Static factory for instantiation of algtool objects
* static ToolFactory<MyRecon> afactory;
* const IToolFactory& RegisterMyRecon = afactory;
* 
* // Standard Constructor
* MyRecon::MyRecon(const std::string& type, 
*                          const std::string& name, 
*                          const IInterface* parent)
*                          : AlgTool( type, name, parent ) {
*     
*     // Declare additional interface
*     declareInterface<IReconTool>(this);
*         
* }
* 
* 
* virtual StatusCode MyRecon::get(const std::string& name, double& x)
* {
*     
*     return StatusCode::SUCCESS; 
* } 
*</pre>
*/


class   IReconTool : virtual public IAlgTool {
public:
 

    // This returns a quantity of interest for reconstruction
    virtual StatusCode get(const std::string& name, double&) = 0;


    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IReconTool; }

};



#endif  // _H_IGlastDetSvc_
