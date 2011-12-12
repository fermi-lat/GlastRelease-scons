/** @file test_EventIntegrityAlg.cxx
    @brief declartion, implementaion of the class test_EventIntegrityAlg

    $Header$
*/
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"


#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "LdfEvent/EventSummaryData.h"

#include <cassert>

/** @class test_EventIntegrityAlg
@brief A simple algorithm.
*/
class test_EventIntegrityAlg : public Algorithm {
public:
    test_EventIntegrityAlg(const std::string& name, ISvcLocator* pSvcLocator);
    /// set parameters and attach to various perhaps useful services.
    StatusCode initialize();
    /// process one event
    StatusCode execute();
    /// clean up
    StatusCode finalize();
    
private: 
    /// number of times called
    int m_count; 
};


//static const AlgFactory<test_EventIntegrityAlg>  Factory;
//const IAlgFactory& test_EventIntegrityAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(test_EventIntegrityAlg);

test_EventIntegrityAlg::test_EventIntegrityAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
    
}

StatusCode test_EventIntegrityAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    return sc;
}

StatusCode test_EventIntegrityAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    m_count++;
    // Retrieving pointers from the TDS
    
    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);

    LdfEvent::EventSummaryData *summary = new LdfEvent::EventSummaryData(0);
    if (m_count % 2) {
        summary->initEventFlags(1);
    } else {
        summary->initEventFlags(0);
    }

    sc = eventSvc()->registerObject("/Event/EventSummary", summary);


    return sc;
}

StatusCode test_EventIntegrityAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    return sc;
}

