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
#include "Event/Digi/AcdDigi.h"
#include "LdfEvent/DiagnosticData.h"
#include "LdfEvent/EventSummaryData.h"
#include "LdfEvent/Gem.h"
#include "LdfEvent/LsfMetaEvent.h"
#include "AdfEvent/AdfEvent.h"
#include "idents/CalXtalId.h"

//static const AlgFactory<test_LdfConverterAlg> Factory;
//const IAlgFactory& test_LdfConverterAlgFactory = Factory

//#include <map>

/** @class test_LdfConverterAlg
 * @brief Takes data from the TDS to test reading from EBF files
 *
 * @author Heather Kelly
 * $Header$
 */

class test_LdfConverterAlg : public Algorithm
{	
public:
    
    test_LdfConverterAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:

    //StatusCode readDigiData();
    StatusCode readEventSummaryData();
    StatusCode readTkrDigiData();
    StatusCode readCalDigiData();
    StatusCode readAcdDigiData();
    StatusCode readAdfData() ;

    int m_count;

    IntegerProperty m_header, m_metaEvent, m_summary, m_gem;
    IntegerProperty m_tkr, m_cal, m_acd, m_diagnostic, m_adf;
    
};

//static const AlgFactory<test_LdfConverterAlg>  Factory;
//const IAlgFactory& test_LdfConverterAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(test_LdfConverterAlg);


test_LdfConverterAlg::test_LdfConverterAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    declareProperty("PrintEventHeader", m_header = 1);
    declareProperty("PrintMetaEvent", m_metaEvent = 1);
    declareProperty("PrintEventSummary", m_summary = 0);
    declareProperty("PrintGem", m_gem = 1);
    declareProperty("PrintTkr", m_tkr = 1);
    declareProperty("PrintCal", m_cal = 1);
    declareProperty("PrintAcd", m_acd = 1);
    declareProperty("PrintAdf", m_adf = 0);
    declareProperty("PrintDiagnostic", m_diagnostic = 0);
}

StatusCode test_LdfConverterAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    m_count = 0;
    
    setProperties();
    
    return sc;
    
}

StatusCode test_LdfConverterAlg::execute()
{

    ++m_count;

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

   if (m_metaEvent) {
       SmartDataPtr<LsfEvent::MetaEvent> meta(eventSvc(), "/Event/MetaEvent");
       if (!meta) 
           log << MSG::INFO << "No MetaEvent on TDS" << endreq;
       else {
           log << MSG::INFO;
           meta->fillStream(log.stream());
           log << endreq;
       }
   }


   if (m_header) {
      SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
      if (!evt) {
         log << MSG::ERROR << "Did not retrieve event" << endreq;
          return StatusCode::FAILURE;
      }
      int eventId = evt->event();
      int run = evt->run();
      log << MSG::INFO << "runId: " << run << " eventId: " << eventId << endreq;
      double time = evt->time();
      log << MSG::INFO << "Time (Seconds since Mission Start): " << time 
          << endreq;
      log << MSG::INFO << "Livetime (seconds) " << evt->livetime() << endreq;

   }

   if (m_gem) {
      SmartDataPtr<LdfEvent::Gem> gem(eventSvc(), "/Event/Gem");
      if (!gem) {
        log << MSG::INFO << "No GEM available" << endreq;
      } else {
        log << MSG::INFO;
        (*gem).fillStream(log.stream());
        log << endreq;
      }
   }
  

    if (m_summary) {
        sc = readEventSummaryData();
        if (sc.isFailure()) {
            return sc;
        }
    }

    if (m_tkr) {
        sc = readTkrDigiData();
        if (sc.isFailure()) 
            return sc;
    }

    if (m_cal) {
        sc = readCalDigiData();
        if (sc.isFailure()) 
            return sc;
    }

    if (m_acd) {
        sc = readAcdDigiData();
        if (sc.isFailure()) 
            return sc;
    }


    if (m_adf) {
        sc = readAdfData();

        if (sc.isFailure()) 
            return sc;
    }

   if (m_diagnostic) {
        SmartDataPtr<LdfEvent::DiagnosticData> diag(eventSvc(), 
                                                    "/Event/Diagnostic");
        if (!diag) {
            log << MSG::ERROR << "Did not retrieve diagnostic data" << endreq;
            return StatusCode::FAILURE;
        }
        // Check diagnostic Data
        log << MSG::INFO;
        diag->fillStream(log.stream());
        log << endreq;
    }

    return sc;
}

StatusCode test_LdfConverterAlg::readEventSummaryData() {
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

StatusCode test_LdfConverterAlg::readTkrDigiData() {
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

StatusCode test_LdfConverterAlg::readCalDigiData() {
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

StatusCode test_LdfConverterAlg::readAcdDigiData() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::DigiEvent> digiEvt(eventSvc(), EventModel::Digi::Event);
    if (!digiEvt) {
        log << MSG::ERROR << "Did not retrieve DigiEvent" << endreq;
        return StatusCode::FAILURE;
    } else {
        log << MSG::INFO << "DigiEvent found!" << endreq;
    }
    
    SmartDataPtr<Event::AcdDigiCol> digiCol(eventSvc(), EventModel::Digi::AcdDigiCol);
    if (digiCol == 0) {
        log << "no AcdDigiCol found" << endreq;
        return sc;
    } else {
        log << " Event " << m_count << ", " << digiCol->size() << " ACD digis found " << endreq;
        int ndigi = 0;

        Event::AcdDigiCol::const_iterator pAcdDigi = digiCol->begin();
        for (pAcdDigi; pAcdDigi!= digiCol->end(); pAcdDigi++) {
            log <<MSG::INFO << "Digi " << ndigi++ << " ";
            (**pAcdDigi).fillStream(log.stream());
            log << endreq;
        }


    }
    return StatusCode::SUCCESS;
}

StatusCode test_LdfConverterAlg::readAdfData() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<DataObject> ancEvt(eventSvc(), "/Event/AncillaryEvent");
    if (!ancEvt) {
        log << MSG::ERROR << "Did not retrieve AncillaryEvent" << endreq;
        return StatusCode::FAILURE;
    }
    SmartDataPtr<AncillaryData::AdfEvent> adfEvt(eventSvc(), "/Event/AncillaryEvent/AdfEvent");
    if (!adfEvt) {
        log << MSG::ERROR << "Did not retrieve AdfEvent" << endreq;
        return StatusCode::FAILURE;
    } else {
        log << MSG::INFO << "AdfEvent found!" << endreq;
    }
    adfEvt->print();
    adfEvt->dump();
    return StatusCode::SUCCESS;
}


StatusCode test_LdfConverterAlg::finalize()
{    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}

