/** @file InterleaveAlg.cxx

@brief declaration and definition of the class InterleaveAlg

$Header$

*/


#include "InterleaveAlg.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"

// access to the tuple
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include "CLHEP/Random/RandFlat.h"

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


double InterleaveAlg::s_triggerRate=1;

//------------------------------------------------------------------------
//! ctor
InterleaveAlg::InterleaveAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
, m_selector(0)
, m_count(0), m_downlink(0),  m_meritTuple(0)
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

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get a pointer to RootTupleSvc
    sc = service("RootTupleSvc", m_rootTupleSvc);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "failed to get the RootTupleSvc" << endreq;
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
    m_magLatLeaf = m_meritTuple->GetLeaf("PtMagLat");
    if(m_magLatLeaf==0) {
        log << MSG::ERROR << "PtMagLat leaf not found in the tuple" << endreq;
        return StatusCode::FAILURE;
    }

    // initialize the background selection
    try {
        std::string file(m_rootFile.value());
        facilities::Util::expandEnvVar(&file);
        m_selector = new BackgroundSelection(file, m_treeName, m_meritTuple);

    }catch( const std::exception& e){
        log << MSG::ERROR << e.what() << endl;
        return StatusCode::FAILURE;
    }

    // these will not be copied from the old tuple.
    m_selector->disable("EvtRun");       
    m_selector->disable("EvtEventId");   
    m_selector->disable("EvtElapsedTime");
    m_selector->disable("EvtLiveTime");  
    m_selector->disable("Pt*");


    // set initial default values
    s_triggerRate = m_selector->triggerRate(0.);
    log << MSG::INFO << "initialized OK: initial trigger rate is " << s_triggerRate << " Hz"<< endreq;

    return sc;
}


//------------------------------------------------------------------------
//! process an event
StatusCode InterleaveAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    setFilterPassed(false); // since this is on a branch, and we want the sequence to fail
    // check that the TDS has an appropriate pseudo-background 

    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);

    if (0==particles) {
        log << MSG::ERROR << "No MC particles!" << endreq;
        return StatusCode::FAILURE;
    }   

    const Event::McParticle& primary = **particles->begin();
    double ke = primary.initialFourMomentum().e()-primary.initialFourMomentum().m();
    log << MSG::DEBUG << "Primary particle energy: " << ke << endreq;

    if( ke>1. ) return sc; // not a flagged sampled_background 
    ++m_count;

    // Found a sampled_background particle
    // perhaps modify its rate
    s_triggerRate = m_selector->triggerRate(magneticLatitude());
    double fraction = m_selector->downlinkRate(magneticLatitude())/s_triggerRate;
    double r = RandFlat::shoot();
    if( r>fraction ) return sc;

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
    static bool done = false;
    if( done ) return sc;
    done=true;
    MsgStream log(msgSvc(), name());

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

