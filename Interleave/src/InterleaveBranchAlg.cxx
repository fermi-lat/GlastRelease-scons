/** @file InterleaveBranchAlg.cxx

@brief declaration and definition of the class InterleaveBranchAlg

$Header$

*/
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
//#include "Event/TopLevel/MCEvent.h"

// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

#include "facilities/Util.h"

/** @class InterleaveBranchAlg
    @brief Gaudi Alg whose sole purpose is to direct interleave events to the interleave sequencer branch
 */


class InterleaveBranchAlg : public Algorithm {
public:
    InterleaveBranchAlg(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    int m_count;       ///< number of processed events
};

//------------------------------------------------------------------------

//static const AlgFactory<InterleaveBranchAlg>  Factory;
//const IAlgFactory& InterleaveBranchAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(InterleaveBranchAlg);

//------------------------------------------------------------------------
//! ctor
InterleaveBranchAlg::InterleaveBranchAlg(const std::string& name, ISvcLocator* pSvcLocator)
              :Algorithm(name, pSvcLocator), m_count(0)
{
    return;
}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode InterleaveBranchAlg::initialize()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode InterleaveBranchAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream   log( msgSvc(), name() );

    // check that the TDS has an appropriate pseudo-background 

    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);

    if (0==particles) 
    {
        log << MSG::ERROR << "No MC particles!" << endreq;
        return StatusCode::FAILURE;
    }   

    const Event::McParticle& primary = **particles->begin();
    double ke = primary.initialFourMomentum().e()-primary.initialFourMomentum().m();

    if( ke > 0. )
    {
        return sc; // not a flagged sampled_background 
    }
    ++m_count;

    // We want to read in a background event so "fail" the filter to direct to the interleave branch
    setFilterPassed(false); 

    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode InterleaveBranchAlg::finalize()
{
    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "Directed "<< m_count << " events to the interleave branch." << endreq; 

    return sc;
}

