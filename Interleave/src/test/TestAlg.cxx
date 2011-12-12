/** @file TestAlg.cxx
    @brief Used for test program

 $Header$

*/

// Include files
#include "FluxSvc/IFluxSvc.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

// access to the tuple
#include "ntupleWriterSvc/INTupleWriterSvc.h"

//#include "../BackgroundSelection.h"
// Gaudi system includes
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"


#include <string>

/*! \class TestAlg
\brief 

*/

class TestAlg : public Algorithm {

public:
    //! Constructor of this form must be provided
    TestAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();


private:
    IFluxSvc* m_fsvc; /// pointer to the flux Service 
    

    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc * m_rootTupleSvc;

    // event id stuff to set, since normally done in AnalysisNtuple.
    unsigned int   EvtRun;
    unsigned int EvtEventId;
    unsigned long long EvtEventId64;
    double EvtElapsedTime;
    float EvtLiveTime;

    float McSourceId; // for testing

    float EvtEnergyCorr
        , TkrNumTracks
        , VtxXDir, VtxYDir, VtxZDir
        , VtxX0, VtxY0, VtxZ0
        , Tkr1XDir, Tkr1YDir, Tkr1ZDir
        , Tkr1X0, Tkr1Y0, Tkr1Z0
        , Tkr1FirstLayer
        , CTBBestEnergy
        , CTBBestXDir, CTBBestYDir, CTBBestZDir
        , CTBBestEnergyProb, CTBBestEnergyRatio, CTBCORE
        , CTBClassLevel, CTBParticleType
        , Cal1MomXDir, Cal1MomYDir, Cal1MomZDir
        ;


};


//static const AlgFactory<TestAlg>  Factory;
//const IAlgFactory& CRTestAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(TestAlg);

//------------------------------------------------------------------------------
//
TestAlg::TestAlg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator)
, EvtRun(99)
, EvtEventId(0)
, EvtEventId64(0)
, EvtElapsedTime(0)
, EvtLiveTime(0)
, McSourceId(0)
, TkrNumTracks(1.)
, CTBBestEnergy(1000.)
, Tkr1FirstLayer (6)
, CTBBestXDir(0), CTBBestYDir(0), CTBBestZDir(1.)
, Cal1MomXDir(0.), Cal1MomYDir(0.), Cal1MomZDir(1.)
{

}

//------------------------------------------------------------------------------
/*! */
StatusCode TestAlg::initialize() {


    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get the service
    StatusCode sc = service("FluxSvc", m_fsvc);
    m_fsvc->GPSinstance()->notifyObservers();

    // get a pointer to RootTupleSvc
    sc = service("RootTupleSvc", m_rootTupleSvc);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "failed to get the RootTupleSvc" << endreq;
        return sc;
    }

    std::string treename("MeritTuple");

    // a minimal set of items to test mechanism.
    m_rootTupleSvc->addItem(treename, "EvtRun",         &EvtRun );
    m_rootTupleSvc->addItem(treename, "EvtEventId",     &EvtEventId );
    m_rootTupleSvc->addItem(treename, "EvtEventId64",   &EvtEventId64);
    m_rootTupleSvc->addItem(treename, "EvtElapsedTime", &EvtElapsedTime );
    m_rootTupleSvc->addItem(treename, "EvtLiveTime",    &EvtLiveTime );
    m_rootTupleSvc->addItem(treename, "McSourceId",     &McSourceId );
    m_rootTupleSvc->addItem(treename, "EvtEnergyCorr" , &EvtEnergyCorr );
    m_rootTupleSvc->addItem(treename, "TkrNumTracks",   &TkrNumTracks);

    m_rootTupleSvc->addItem(treename,"VtxXDir",            &VtxXDir     );
    m_rootTupleSvc->addItem(treename,"VtxYDir",            &VtxYDir     );    
    m_rootTupleSvc->addItem(treename,"VtxZDir",            &VtxZDir     );
    m_rootTupleSvc->addItem(treename,"VtxX0",              &VtxX0       );
    m_rootTupleSvc->addItem(treename,"VtxY0",              &VtxY0       );
    m_rootTupleSvc->addItem(treename,"VtxZ0",              &VtxZ0       );
    m_rootTupleSvc->addItem(treename,"Tkr1XDir",           &Tkr1XDir    );
    m_rootTupleSvc->addItem(treename,"Tkr1YDir",           &Tkr1YDir    );
    m_rootTupleSvc->addItem(treename,"Tkr1ZDir",           &Tkr1ZDir    );
    m_rootTupleSvc->addItem(treename,"Tkr1X0",             &Tkr1X0      );
    m_rootTupleSvc->addItem(treename,"Tkr1Y0",             &Tkr1Y0      );
    m_rootTupleSvc->addItem(treename,"Tkr1Z0",             &Tkr1Z0      );
    m_rootTupleSvc->addItem(treename,"Tkr1FirstLayer",     &Tkr1FirstLayer);
    m_rootTupleSvc->addItem(treename,"CTBBestEnergy",      &CTBBestEnergy);
    m_rootTupleSvc->addItem(treename,"CTBBestXDir",        &CTBBestXDir );  
    m_rootTupleSvc->addItem(treename,"CTBBestYDir",        &CTBBestYDir );  
    m_rootTupleSvc->addItem(treename,"CTBBestZDir",        &CTBBestZDir );  
    m_rootTupleSvc->addItem(treename,"CTBBestEnergyProb",  &CTBBestEnergyProb );  
    m_rootTupleSvc->addItem(treename,"CTBBestEnergyRatio", &CTBBestEnergyRatio );
    m_rootTupleSvc->addItem(treename,"CTBCORE",            &CTBCORE );
    m_rootTupleSvc->addItem(treename,"CTBClassLevel",      &CTBClassLevel );
    m_rootTupleSvc->addItem(treename,"CTBParticleType",    &CTBParticleType );
    m_rootTupleSvc->addItem(treename,"Cal1MomXDir",        &Cal1MomXDir );
    m_rootTupleSvc->addItem(treename,"Cal1MomYDir",        &Cal1MomYDir );
    m_rootTupleSvc->addItem(treename,"Cal1MomZDir",        &Cal1MomZDir );

    return sc;
}


//------------------------------------------------------------------------------
StatusCode TestAlg::execute() {

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    

    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);

    EvtRun           = header->run();
    EvtEventId       = header->event();
    EvtElapsedTime   = header->time();
    EvtLiveTime      = header->livetime();
    
    return sc;
}


//------------------------------------------------------------------------------
StatusCode TestAlg::finalize() {

    return StatusCode::SUCCESS;
}






