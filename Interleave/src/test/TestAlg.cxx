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

#include "../BackgroundSelection.h"
// Gaudi system includes
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"


#include <string>
// copy here for testing.
//------------------------------------------------------------------------
double BackgroundSelection::triggerRate(double maglat)
{
    //y = 0.0481x3 - 0.9295x2 - 27.72x + 2082.2
    static double c[]={2082.2, -27.72, -0.9295, 0.0481};

    double x = fabs(maglat);
    if( x< 8) x=8;
    return  ((c[3]*x + c[2])*x +c[1])*x + c[0];
}

//------------------------------------------------------------------------
double BackgroundSelection::downlinkRate(double  maglat )
{
    //y = 0.0094x3 - 0.1951x2 - 5.4218x + 508.48
    static double c[]={508.48, -5.4218, -0.1951, 0.0094};

    double x = fabs(maglat);
    if( x< 8) x=8;
    return ((c[3]*x + c[2])*x +c[1])*x + c[0];

}
//----

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
    float EvtRun;
    float EvtEventId;
    double EvtElapsedTime;
    float EvtLiveTime;

    float McSourceId; // for testing

};


static const AlgFactory<TestAlg>  Factory;
const IAlgFactory& CRTestAlgFactory = Factory;

//------------------------------------------------------------------------------
//
TestAlg::TestAlg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator)
, EvtRun(99)
, EvtEventId(0)
, EvtElapsedTime(0)
, EvtLiveTime(0)
, McSourceId(0)
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
    m_rootTupleSvc->addItem(treename, "EvtElapsedTime", &EvtElapsedTime );
    m_rootTupleSvc->addItem(treename, "EvtLiveTime",    &EvtLiveTime );
    m_rootTupleSvc->addItem(treename, "McSourceId",     &McSourceId );

    std::cout << "Testing the rates\n\t lat     trigger     downlink\n"<< std::endl;
    for( double lat = 0; lat<45; lat+=5){
        std::cout << "\t" << lat << " "<< BackgroundSelection::triggerRate(lat) << "   " << BackgroundSelection::downlinkRate(lat) << std::endl;
    }
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






