#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/TkrDigi.h"
#include "Event/Digi/CalDigi.h"
#include "LdfEvent/DiagnosticData.h"
#include "LdfEvent/EventSummaryData.h"
#include "idents/CalXtalId.h"


//#include <map>

/** @class test_EbfConverterAlg
 * @brief Takes data from the TDS to test reading from EBF files
 *
 * @author Heather Kelly
 * $Header$
 */

class test_EbfConverterAlg : public Algorithm
{	
public:
    
    test_EbfConverterAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:

    //StatusCode readDigiData();
    StatusCode readEventSummaryData();
    StatusCode readTkrDigiData();
    StatusCode readCalDigiData();

    int m_count;
    
};

static const AlgFactory<test_EbfConverterAlg>  Factory;
const IAlgFactory& test_EbfConverterAlgFactory = Factory;


test_EbfConverterAlg::test_EbfConverterAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
}

StatusCode test_EbfConverterAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    m_count = 0;
    
    setProperties();
    
    return sc;
    
}

StatusCode test_EbfConverterAlg::execute()
{

    ++m_count;

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

  //  sc = readDigiData();

    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
    if (!evt) {
        log << MSG::ERROR << "Did not retrieve event" << endreq;
        return StatusCode::FAILURE;
    }
    //Check event header
    int eventId = evt->event();
    int run = evt->run();
    log << MSG::INFO << "runId: " << run << " eventId: " << eventId << endreq;

    /*SmartDataPtr<LdfEvent::DiagnosticData> diag(eventSvc(), "/Event/Diagnostic");
    if (!diag) {
        log << MSG::ERROR << "Did not retrieve diagnostic data" << endreq;
        return StatusCode::FAILURE;
    }
    // Check diagnostic Data
    log << MSG::INFO;
    diag->fillStream(log.stream());
    log << endreq;
*/

    sc = readEventSummaryData();
    if (sc.isFailure()) {
        return sc;
    }

    sc = readTkrDigiData();
    if (sc.isFailure()) {
        return sc;
    }

    sc = readCalDigiData();

    return sc;
}

StatusCode test_EbfConverterAlg::readEventSummaryData() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<LdfEvent::EventSummaryData> summary(eventSvc(), "/Event/EventSummary");
    if (!summary) {
       log << MSG::ERROR << "Did not retrieve event summary" << endreq;
       return StatusCode::FAILURE;
    } else {
       log << MSG::INFO << "Event Summary found!";
       (*summary).fillStream(log.stream());
       log << endreq;
    }

    return sc;
}

StatusCode test_EbfConverterAlg::readTkrDigiData() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::DigiEvent> digiEvt(eventSvc(), EventModel::Digi::Event);
    if (!digiEvt) {
        log << MSG::ERROR << "Did not retrieve DigiEvent" << endreq;
        return StatusCode::FAILURE;
    } else {
        log << MSG::INFO << "DigiEvent found!" << endreq;
    }
    
    SmartDataPtr<Event::TkrDigiCol> digiCol(eventSvc(), EventModel::Digi::TkrDigiCol);

    if (digiCol == 0) {
        log << "no TkrDigiCol found" << endreq;
        return sc;
    } else {
        log << " Event " << m_count << ", " << digiCol->size() << " TKR digis found " << endreq;
        if(m_count<10) {
            log << MSG::INFO << endreq << "Detailed dump of 1st event: " << endreq << endreq;
            int ndigi = 0;
            Event::TkrDigiCol::const_iterator pTkrDigi = digiCol->begin();
            for (pTkrDigi; pTkrDigi!= digiCol->end(); pTkrDigi++) {
                log <<MSG::INFO << "Digi " << ndigi++ << " ";
                // output the digi proper:  " << digi " no longer works because
                // "digi" now returns an integer sort order
                (**pTkrDigi).fillStream(log.stream());
                log << endreq;
            }
        }
    }
    

    return StatusCode::SUCCESS;
}

StatusCode test_EbfConverterAlg::readCalDigiData() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::DigiEvent> digiEvt(eventSvc(), EventModel::Digi::Event);
    if (!digiEvt) {
        log << MSG::ERROR << "Did not retrieve DigiEvent" << endreq;
        return StatusCode::FAILURE;
    } else {
        log << MSG::INFO << "DigiEvent found!" << endreq;
    }
    
    SmartDataPtr<Event::CalDigiCol> digiCol(eventSvc(), EventModel::Digi::CalDigiCol);

    if (digiCol == 0) {
        log << "no CalDigiCol found" << endreq;
        return sc;
    } else {
        log << " Event " << m_count << ", " << digiCol->size() << " CAL digis found " << endreq;
        int ndigi = 0;

        Event::CalDigiCol::const_iterator pCalDigi = digiCol->begin();
        for (pCalDigi; pCalDigi!= digiCol->end(); pCalDigi++) {
            log <<MSG::INFO << "Digi " << ndigi++ << " ";
            // output the digi proper:  " << digi " no longer works because
            // "digi" now returns an integer sort order
            (**pCalDigi).fillStream(log.stream());
            log << endreq;
        }

    }

    return StatusCode::SUCCESS;

}

/*
StatusCode test_EbfConverterAlg::readDigiData() {
MsgStream log(msgSvc(), name());
StatusCode sc = StatusCode::SUCCESS;

SmartDataPtr<Event::CalDigiCol> calDigiColTds(eventSvc(), EventModel::Digi::CalDigiCol);
if (!calDigiColTds) return sc;
Event::CalDigiCol::const_iterator calDigiTds;

for (calDigiTds = calDigiColTds->begin(); calDigiTds != calDigiColTds->end(); calDigiTds++) {
log << MSG::DEBUG;
(*calDigiTds)->fillStream(log.stream());
log << endreq;
idents::CalXtalId::CalTrigMode modeTds = (*calDigiTds)->getMode();
idents::CalXtalId idTds = (*calDigiTds)->getPackedId();
if (modeTds == idents::CalXtalId::BESTRANGE) {
const Event::CalDigi::CalXtalReadout *readoutTds = (*calDigiTds)->getXtalReadout(0);
} else {
int range;
for (range = idents::CalXtalId::LEX8; range <= idents::CalXtalId::HEX1; range++) {
const Event::CalDigi::CalXtalReadout *readoutTds = (*calDigiTds)->getXtalReadout(range);
}
}
}

return sc;
}
*/


StatusCode test_EbfConverterAlg::finalize()
{    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}

