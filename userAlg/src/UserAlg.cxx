// $Header$

// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

// ntuple interface
#include "ntupleWriterSvc/INTupleWriterSvc.h"

// if using the gui
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
#include "gui/DisplayControl.h"

// TDS class declarations: input data, and McParticle tree
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/MonteCarlo/McIntegratingHit.h"
#include "GlastEvent/MonteCarlo/McPositionHit.h"

#include "GlastEvent/TopLevel/Event.h"
#include "GlastEvent/TopLevel/EventModel.h"

// for access to instrument.ini
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "xml/IFile.h"

#include <cassert>

/** @class UserAlg
@brief A simple algorithm.
*/
class UserAlg : public Algorithm {
public:
    UserAlg(const std::string& name, ISvcLocator* pSvcLocator);
    /// set parameters and attach to various perhaps useful services.
    StatusCode initialize();
    /// process one event
    StatusCode execute();
    /// clean up
    StatusCode finalize();
    
private: 
    /// number of times called
    int m_count; 
    /// the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc; 
    /// access to the Gui Service for display of 3-d objects
    IGuiSvc*    m_guiSvc;
    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc *m_ntupleWriteSvc;
    /// parameter to store the logical name of the ROOT file to write to
    std::string m_tupleName;
};


static const AlgFactory<UserAlg>  Factory;
const IAlgFactory& UserAlgFactory = Factory;

UserAlg::UserAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0), m_detSvc(0), m_guiSvc(0)
{
    // declare properties with setProperties calls
    declareProperty("tupleName",  m_tupleName="");
    
}

StatusCode UserAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    if( m_tupleName.empty()) {
        log << MSG::INFO << "tupleName property not set!  No ntuple output"<<endreq;
    }
    
    // now try to find the GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc).isFailure()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    
    // get the Gui service (not required)
    if (service("GuiSvc", m_guiSvc).isFailure ()){
        log << MSG::WARNING << "No GuiSvc, so no display" << endreq;
    }else{
        //m_guiSvc->guiMgr()->display().add(new Rep, "User rep");
    }
    
    // get a pointer to our ntupleWriterSvc
    if (!m_tupleName.empty()) {
        if (service("ntupleWriterSvc", m_ntupleWriteSvc, true).isFailure()) {
            log << MSG::ERROR << "UserAlg failed to get the ntupleWriterSvc" << endreq;
            return StatusCode::FAILURE;
        }
    }
    return sc;
}

StatusCode UserAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "executing " << ++m_count << " time" << endreq;
    
    // Here we are adding to our ROOT ntuple
    if (!m_tupleName.empty())
        sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "NumCalls", m_count);
    
    // An example of retrieving data from the TDS
    SmartDataPtr<mc::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);
    
    if (particles) {
        
        mc::McParticleCol::const_iterator p;
        
        for (p = particles->begin(); p != particles->end(); p++) {
            
        }
    }
    


    // An example of pausing the display 
    bool pause = false;
    // calculate some condition here...
    // pause = .....;
    if (pause) {
        log << MSG::INFO << "Pausing at event " << m_count <<  endreq;
        m_guiSvc->guiMgr()->pause();
    } 
    
    
    return sc;
}

StatusCode UserAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
    return sc;
}

