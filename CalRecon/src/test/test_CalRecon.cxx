
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

class test_CalRecon : public Algorithm {
public:
    test_CalRecon(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    //! number of times called
    int m_count;
    
};
//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_CalRecon );

static const AlgFactory<test_CalRecon>  Factory;
const IAlgFactory& test_CalReconFactory = Factory;

//------------------------------------------------------------------------
//! ctor
test_CalRecon::test_CalRecon(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode test_CalRecon::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_CalRecon::execute() {
    
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "executing " << ++m_count << " time" << endreq;


    // First, the collection of CalRecons is retrieved from the TDS
    SmartDataPtr<Event::CalClusterCol> clusCol(eventSvc(),EventModel::CalRecon::CalClusterCol) ;

    if ((clusCol == 0) || (clusCol->size()) == 0) {
        log << MSG::FATAL << "no/incorrect number calorimeter clusters found" << endreq;
        return StatusCode::FAILURE ;
    }
    
    // David: I display here what is not already displayed
    // by the CalClustersAlg
    Event::CalClusterCol::const_iterator cluster ;
    for ( cluster = clusCol->begin() ;
          cluster != clusCol->end() ;
          ++cluster )
     {    
      (*cluster)->writeOut(log<<MSG::INFO) ; 
      log<<MSG::INFO<<"Profile Corr Energy "<<(*cluster)->getFitEnergy()<<endreq ; 
      log<<MSG::INFO<<"Last Layer Corr Energy "<<(*cluster)->getEnergyLeak()<<endreq ; 
     }

    return StatusCode::SUCCESS ;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_CalRecon::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;
    
    return sc;
}


