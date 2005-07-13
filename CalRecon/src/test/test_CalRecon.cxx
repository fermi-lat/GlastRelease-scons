
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include <cmath>
#include <cerrno>
#include <exception>

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
    MsgStream log(msgSvc(),name()) ;
    log<<MSG::INFO<<"initialize"<<endreq ;
    
    // check errno
    errno = 0 ;
    log10(-1.) ;
    if (!errno) {
        log<<MSG::FATAL<<"errno not workng as expected"<<endreq ;
        return StatusCode::FAILURE ;
    }
    
    return StatusCode::SUCCESS ;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_CalRecon::execute() {
    
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "executing " << ++m_count << " time" << endreq;


    // First, the collection of CalRecons is retrieved from the TDS
    SmartDataPtr<Event::CalClusterCol> clusCol(eventSvc(),EventModel::CalRecon::CalClusterCol) ;
    SmartDataPtr<Event::CalEventEnergy> calEnergy(eventSvc(),EventModel::CalRecon::CalEventEnergy);

    if ((clusCol == 0) || (clusCol->size()) == 0) {
        throw std::runtime_error("no/incorrect number calorimeter clusters found") ;
    }
    
    Event::CalClusterCol::const_iterator cluster ;
    for ( cluster = clusCol->begin() ;
          cluster != clusCol->end() ;
          ++cluster )
    {
      log<<MSG::INFO<<"Energy "<<(*cluster)->getCalParams().getEnergy()<<endreq ;
      log<<MSG::INFO<<"NbTruncXtals " <<(*cluster)->getNumTruncXtals()<<endreq ;
      log<<MSG::INFO<<"Position"
        <<" "<<(*cluster)->getPosition().x() 
        <<" "<<(*cluster)->getPosition().y() 
        <<" "<<(*cluster)->getPosition().z()
        <<endreq ;
      log<<MSG::INFO<<"Direction"
        <<" "<<(*cluster)->getDirection().x() 
        <<" "<<(*cluster)->getDirection().y() 
        <<" "<<(*cluster)->getDirection().z()
        <<endreq ;
    }

    Event::CalCorToolResultCol::iterator corIter ;
    for( corIter= calEnergy->begin(); corIter != calEnergy->end(); ++corIter )
    {
        Event::CalCorToolResult * corResult = *corIter;
        if (corResult->getCorrectionName() == "CalProfileTool" )
        {
            log<<MSG::INFO<<"Profile Corr Energy "
              <<corResult->getParams().getEnergy()
              <<endreq ; 
        }
        else if (corResult->getCorrectionName() == "CalFullProfileTool" )
        {
            log<<MSG::INFO<<"FullProfile Corr Energy "
	       <<corResult->getParams().getEnergy()
	       << " alpha = " << (*corResult)["alpha"] 
	       << " tmax = " << (*corResult)["tmax"] 
	       << " tkrrln = " << (*corResult)["tkrrln"] 
	       << " lastx0 = " << (*corResult)["lastx0"] 
	       << " totchisq = " << (*corResult)["totchisq"] 
	       << " chisq = " << (*corResult)["chisq"] 
	       << " parcf = " << (*corResult)["parcf"] 
	       << " parc = " << (*corResult)["parc"] 
              <<endreq ; 
        }
        else if (corResult->getCorrectionName() == "CalLastLayerLikelihoodTool" )
        {
            double llStatus = (*corResult)["llStatus"] ;
            if (llStatus==0)
              log<<MSG::INFO<<"Last Layer Corr Energy "
                <<corResult->getParams().getEnergy()
                <<endreq ; 
            else
              log<<MSG::INFO<<"Last Layer Corr Energy "
                <<llStatus
                <<endreq ; 
        }
        else if (corResult->getCorrectionName() == "CalValsCorrTool" )
        {
            log<<MSG::INFO<<"CalValsCorrTool Energy "
              <<corResult->getParams().getEnergy()
              <<endreq ; 
        }
        else
        {
            log<<MSG::WARNING<<"UNKNOWN CORRECTION NAME: "
              <<corResult->getCorrectionName()
              <<endreq ; 
        }
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


