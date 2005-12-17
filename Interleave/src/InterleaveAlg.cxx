/** @file InterleaveAlgcxx

@brief declaration and definition of the class InterleaveAlg

$Header$

*/


#include "InterleaveAlg.h"

// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"

// to write a Tree with entryosure info
#include "ntupleWriterSvc/INTupleWriterSvc.h"
  
// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"


#include <cassert>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

//------------------------------------------------------------------------

static const AlgFactory<InterleaveAlg>  Factory;
const IAlgFactory& InterleaveAlgFactory = Factory;

double InterleaveAlg::s_triggerRate=1;
double InterleaveAlg::s_downlinkRate;

//------------------------------------------------------------------------
//! ctor
InterleaveAlg::InterleaveAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
, m_count(0), m_downlink(0), m_event(0), m_meritTuple(0)
{
    // declare properties with setProperties calls
    declareProperty("TriggerRate",  s_triggerRate=10000.);
    declareProperty("DownLinkRate", s_downlinkRate=400.);
    declareProperty("RootFile",     m_rootFile="");
    declareProperty("TreeName",     m_treeName="MeritTuple");

}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode InterleaveAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get a pointer to RootTupleSvc
    sc = service("RootTupleSvc", m_rootTupleSvc);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "failed to get the RootTupleSvc" << endreq;
        return sc;
    }
    TFile* m_tfile = new TFile(m_rootFile.value().c_str(), "readonly");
    if( 0==m_tfile || !m_tfile->IsOpen()){
        log << MSG::ERROR << "Could not open the root file " << m_rootFile << endreq;
        return StatusCode::FAILURE;
    }
    m_backgroundTuple =  static_cast<TTree*>(m_tfile->Get(m_treeName.value().c_str()));
    if( 0==m_backgroundTuple) {
        log << MSG::ERROR << "Did not find the tree in the root file " << m_rootFile << endreq;
        return StatusCode::FAILURE;
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
    return sc;
}


//------------------------------------------------------------------------
//! process an event
StatusCode InterleaveAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

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

    // Found a sampled_background particle

    ++m_count;
    ++m_downlink;

    if( ++m_event > m_backgroundTuple->GetEntries() ) {
        log << MSG::WARNING << "Ran out of event to mix: restarting" << endreq;
        m_event = 0;
    }
    m_backgroundTuple->GetEvent(m_event);
    
    // revise stuff in the background tuple to agree with the current event heaader
    copyEventInfo(m_backgroundTuple);


    // now copy every value from the background tree to the local merit tree
    copyTreeData(m_backgroundTuple, m_meritTuple);


    m_rootTupleSvc->storeRowFlag( m_treeName.value(), true);

    return StatusCode::SUCCESS;
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

}
//------------------------------------------------------------------------
void InterleaveAlg::copyTreeData(TTree * in, TTree* out)
{
    // TODO: fill this
    copyEventInfo(out); // the full copy would do this.

}
