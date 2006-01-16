/** @file InterleaveAlg.cxx

@brief declaration and definition of the class InterleaveAlg

$Header$

*/

#include "InterleaveAlg.h"

// this service manages livetime calculations
#include "Trigger/ILivetimeSvc.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"

// access to the tuple
#include "ntupleWriterSvc/INTupleWriterSvc.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"

// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

#include "facilities/Util.h"

#include <cassert>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

//------------------------------------------------------------------------

static const AlgFactory<InterleaveAlg>  Factory;
const IAlgFactory& InterleaveAlgFactory = Factory;


double InterleaveAlg::s_rate=0;

//------------------------------------------------------------------------
//! ctor
InterleaveAlg::InterleaveAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
, m_selector(0)
, m_count(0), m_downlink(0),  m_meritTuple(0)
, m_magLatLeaf(0) 

{
    // declare properties with setProperties calls
    declareProperty("RootFile",     m_rootFile="");
    declareProperty("TreeName",     m_treeName="MeritTuple");

}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode InterleaveAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    if( name()=="interleave_filter") return sc;
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get a pointer to RootTupleSvc
    sc = service("RootTupleSvc", m_rootTupleSvc);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "failed to get the RootTupleSvc" << endreq;
        return sc;
    }

    sc = service("LivetimeSvc", m_LivetimeSvc);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "failed to get the LivetimeSvc" << endreq;
        return sc;
    }

    if( m_meritTuple==0){
        void * ptr= 0;
        m_rootTupleSvc->getItem(m_treeName.value().c_str(),"", ptr);
        if( ptr==0){
            log << MSG::ERROR << "Could not find the merit tuple" << endreq;
            return StatusCode::FAILURE;
        }
        m_meritTuple = static_cast<TTree*>(ptr);
    }
    // initialize the background selection
    try {
        std::string file(m_rootFile.value());
        if( !file.empty()){
            facilities::Util::expandEnvVar(&file);
            log << MSG::INFO << "Using file " << file << " for interleave." << endreq;
            m_selector = new BackgroundSelection(file,  m_meritTuple);
        }else{
            log << MSG::INFO<< "No file specified for interleave" << endreq;
            return sc;
        }

    }catch( const std::exception& e){
        log << MSG::WARNING << e.what() << endreq;
        log << MSG::WARNING << "Continuing without background" << endreq;
            return sc;
        //return StatusCode::FAILURE;
    }

    // these will not be copied from the old tuple.
#if 0 // try it without these two, so that we have access to the original values
    m_selector->disable("EvtRun");       
    m_selector->disable("EvtEventId");   
#endif
    m_selector->disable("EvtElapsedTime");
    m_selector->disable("EvtLiveTime");  
    m_selector->disable("Pt*");
    m_selector->disable("FT1*");
    m_selector->disable("CT*");

    // now modify all leaves to point to our tuple
    m_selector->setLeafPointers();

    // set initial default values for downlink rate to fold in
    s_rate = m_selector->downlinkRate(0.);
    log << MSG::INFO << "initialized OK: initial downlink rate to merge is " << s_rate << " Hz"<< endreq;
    m_LivetimeSvc->setTriggerRate(m_selector->triggerRate(0));

    return sc;
}


//------------------------------------------------------------------------
//! process an event
StatusCode InterleaveAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;

    if( m_selector==0) return sc;
    MsgStream   log( msgSvc(), name() );

    // check that the TDS has an appropriate pseudo-background 

    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);
    
    if( m_magLatLeaf==0){
        if( (m_magLatLeaf = m_meritTuple->GetLeaf("PtMagLat"))==0){
            log << MSG::ERROR << "PtMagLat leaf not found in the tuple" << endreq;
            return StatusCode::FAILURE;
        }
    }

    if (0==particles) {
        log << MSG::ERROR << "No MC particles!" << endreq;
        return StatusCode::FAILURE;
    }   

    const Event::McParticle& primary = **particles->begin();
    double ke = primary.initialFourMomentum().e()-primary.initialFourMomentum().m();
    log << MSG::DEBUG << "Primary particle energy: " << ke << endreq;

    // A special instance of this alg can be used as a filter to kill a sequence
    if( name()=="interleave_filter") {
        if( ke==0) setFilterPassed(false);
        return sc;
    }
    if( ke>1. ){
        setFilterPassed(false); // since this is on a branch, and we want the sequence to fail
        return sc; // not a flagged sampled_background 
    }
    ++m_count;


    m_selector->downlinkRate(magneticLatitude());
    m_LivetimeSvc->setTriggerRate(m_selector->triggerRate(magneticLatitude()));

    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);
    double time  = header->time();

    bool ok = m_LivetimeSvc->isLive(time);

    header->setLivetime( m_LivetimeSvc->livetime());

    ++m_downlink;

    // ask for a tree corresponding to our current position: it will set all the tuple

    m_selector->selectEvent(magneticLatitude());

    // overwrite the event info
    copyEventInfo(m_meritTuple);
    
    // finally flag that we want to add this event to the output tuple
    m_rootTupleSvc->storeRowFlag( m_treeName.value(), true);

    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode InterleaveAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    if( m_selector==0) return sc;
    MsgStream log(msgSvc(), name());

    if( name()=="interleave_filter") {
        return sc;
    }

    log << MSG::INFO << "Processed "<< m_count << " sampled background events, of which "<< m_downlink<< " were passed." << endreq; 
    delete m_selector;
    return sc;
}

//------------------------------------------------------------------------
namespace {
    //! Utility to set a scalar value in a ROOT tree. 
    //! Assume that the type is known!
    template <typename  Type> 
        void setLeafValue(TLeaf* leaf, Type newvalue)
    {
        Type& rval = *static_cast<Type*>(leaf->GetValuePointer());
        rval = newvalue;
    }

}
//------------------------------------------------------------------------
void InterleaveAlg::copyEventInfo(TTree * tree)
{

    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);
    // these types *must* correspond with those in EvtValsTool
    float EvtRun           = header->run();
    float EvtEventId       = header->event();
    double EvtElapsedTime  = header->time();
    float EvtLiveTime      = header->livetime();

    setLeafValue(tree->GetLeaf("EvtRun"),        EvtRun);
    setLeafValue(tree->GetLeaf("EvtEventId"),    EvtEventId);
    setLeafValue(tree->GetLeaf("EvtElapsedTime"),EvtElapsedTime);
    setLeafValue(tree->GetLeaf("EvtLiveTime"),   EvtLiveTime);

    // finally, make the source id negative; make zero -1 by offset. (2's complement)
    float & sourceid = *static_cast<float*>(tree->GetLeaf("McSourceId")->GetValuePointer());
    sourceid=-1-sourceid;

}

//------------------------------------------------------------------------
double InterleaveAlg::magneticLatitude()
{
    return m_magLatLeaf->GetValue();
}

