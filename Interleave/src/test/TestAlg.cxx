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
    float EvtRun;
    float EvtEventId;
    double EvtElapsedTime;
    float EvtLiveTime;

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
{

}

//------------------------------------------------------------------------------
/*! */
StatusCode TestAlg::initialize() {


    MsgStream log(msgSvc(), name());

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

    m_rootTupleSvc->addItem(treename, "EvtRun",         &EvtRun );
    m_rootTupleSvc->addItem(treename, "EvtEventId",     &EvtEventId );
    m_rootTupleSvc->addItem(treename, "EvtElapsedTime", &EvtElapsedTime );
    m_rootTupleSvc->addItem(treename, "EvtLiveTime",    &EvtLiveTime );
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






