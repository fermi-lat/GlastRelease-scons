// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

// ntuple interface
#include "ntupleWriterSvc/INTupleWriterSvc.h"

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
    //! access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc *m_ntupleWriteSvc;
    //! parameter to store the logical name of the ROOT file to write to
    std::string m_tupleName;
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
    declareProperty("tupleName",  m_tupleName="");
    
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

    // get a pointer to our ntupleWriterSvc
    if (service("ntupleWriterSvc", m_ntupleWriteSvc).isFailure()) {
        log << MSG::ERROR << "writeJunkAlg failed to get the ntupleWriterSvc" << endreq;
        return StatusCode::FAILURE;
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

    // Here we are adding to our ROOT ntuple
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "NumCalls", m_count);
 
    // An example of retrieving data from the TDS
    SmartDataPtr<TdGlastData> glastData(eventSvc(),"/Event/TdGlastData");

    // retrieve TKR data pointer
    const SiData *tkrDigiData = glastData->getSiData();
    if (tkrDigiData == 0) {
        log << MSG::INFO << "No TKR Data available" << endreq;
    } else {
        // Get the number of hits in the X and Y layers of plane 0
        int nx = tkrDigiData->nHits(SiData::X, 0);
        int ny = tkrDigiData->nHits(SiData::Y, 0);
        // Now to add the total number of Hits in Plane 0 to the ntuple
        sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "TKR_numHitsPlane0", (nx+ny));
    }

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



