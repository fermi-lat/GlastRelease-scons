// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

// Gaudi ntuple interface
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/INTuple.h"
#include "GaudiKernel/NTuple.h"

// if use the gui
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
#include "gui/DisplayControl.h"

// TDS class declarations: input data, and McParticle tree
#include "GlastEvent/data/TdGlastData.h"
#include "GlastEvent/MonteCarlo/McVertex.h"
#include "GlastEvent/MonteCarlo/McParticle.h"

// for access to instrument.ini
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "xml/IFile.h"

#include <cassert>

// Define the class here instead of in a header file: not needed anywhere but here!
//------------------------------------------------------------------------------
/** 
A simple algorithm.

  
*/
class UserAlg : public Algorithm {
public:
    UserAlg(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    //! number of times called
    int m_count; 
    //! the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc; 
    //! constants from the "instrument.xml" file.
    xml::IFile * m_ini; 
    //! access to the Gui Service for display of 3-d objects
    IGuiSvc*    m_guiSvc;
    
};
//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( UserAlg );

static const AlgFactory<UserAlg>  Factory;
const IAlgFactory& UserAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
UserAlg::UserAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0), m_detSvc(0), m_ini(0), m_guiSvc(0)
{
    // declare properties with setProperties calls
    
}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode UserAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // now try to find the GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc).isFailure()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    // get the ini file
    m_ini = const_cast<xml::IFile*>(m_detSvc->iniFile()); //OOPS!
    assert(4==m_ini->getInt("glast", "xNum"));  // simple check

    // get the Gui service (not required)
    if (service("GuiSvc", m_guiSvc).isFailure ()){
        log << MSG::WARNING << "No GuiSvc, so no display" << endreq;
    }else{
       //m_guiSvc->guiMgr()->display().add(new Rep, "User rep");
    }

    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode UserAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "executing " << ++m_count << " time" << endreq;
    
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode UserAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
    return sc;
}



