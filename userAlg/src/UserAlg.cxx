/** @file UserAlg.cxx
    @brief declartion, implementaion of the class UserAlg

    $Header$
*/
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"


// if using the gui
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
#include "gui/DisplayControl.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"

#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

// for access to geometry perhaps
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

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
    double m_count; 
    /// the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc; 
    /// access to the Gui Service for display of 3-d objects
    IGuiSvc*    m_guiSvc;

    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc *m_rootTupleSvc;

    /// parameter to store the logical name of the ROOT file to write to
    std::string m_treeName;
};


static const AlgFactory<UserAlg>  Factory;
const IAlgFactory& UserAlgFactory = Factory;

UserAlg::UserAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0), m_detSvc(0), m_guiSvc(0)
{
    // declare properties with setProperties calls
    declareProperty("treeName",  m_treeName="");
    
}

StatusCode UserAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    if( m_treeName.empty()) {
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
    
    // example code to create a ROOT tuple and schedule it for filling 

    // get a pointer to RootTupleSvc
    sc = service("RootTupleSvc", m_rootTupleSvc);

    if( sc.isFailure() ) {
        log << MSG::ERROR << "userAlg failed to get the RootTupleSvc" << endreq;
        return sc;
    }
    m_rootTupleSvc->addItem(m_treeName, "count", &m_count);

    return sc;
}

StatusCode UserAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "executing " << ++m_count << " time" << endreq;

    if( m_rootTupleSvc!=0) {
        // necessary to store the tuple entry
        m_rootTupleSvc->storeRowFlag(m_treeName, true);
    }
    

    // Retrieving pointers from the TDS
    
    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);
    SmartDataPtr<Event::MCEvent>     mcheader(eventSvc(),    EventModel::MC::Event);
    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);
    SmartDataPtr<Event::TkrVertexCol>  tracks(eventSvc(),    EventModel::TkrRecon::TkrVertexCol);
    SmartDataPtr<Event::CalClusterCol> clusters(eventSvc(),  EventModel::CalRecon::CalClusterCol);
    SmartDataPtr<Event::CalXtalRecCol> xtalrecs(eventSvc(),  EventModel::CalRecon::CalXtalRecCol);
    SmartDataPtr<Event::AcdRecon> acdrec(eventSvc(),         EventModel::AcdRecon::Event);

    double t = header->time();
    log << MSG::DEBUG << "Event time: " << t << endreq;;

    
    // lop over the monte carlo particle collection
    if (particles) {
        
        Event::McParticleCol::const_iterator piter;
        
        for (piter = particles->begin(); piter != particles->end(); piter++) {
            const Event::McParticle& mcp = **piter;
            
        }
    }
    
#if 0 // enable to use pause with the display
    // An example of pausing the display 
    bool pause = false;
    // calculate some condition here...
    // pause = .....;
    if (pause) {
        log << MSG::INFO << "Pausing at event " << m_count <<  endreq;
        m_guiSvc->guiMgr()->pause();
    } 
#endif
#if 0
    log << MSG::FATAL << "test of a fatal error " << endreq;
    sc = StatusCode::FAILURE;
#endif
    return sc;
}

StatusCode UserAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
    return sc;
}

