/** 
* @file EventIntegrityAlg.cxx
* @brief Declaration and definition of the algorithm EventIntegrityAlg.
*
*  $Header$
*/

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"

#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "LdfEvent/EventSummaryData.h"

/*! \class EventIntegrityAlg
\brief  alg that sets trigger information

@section Attributes for job options:
@param run [0] For setting the run number
@param mask [-1] mask to apply to trigger word. -1 means any, 0 means all.

*/

class EventIntegrityAlg : public Algorithm {

public:

    //! Constructor of this form must be provided
    EventIntegrityAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:

    unsigned int m_mask;

};

static const AlgFactory<EventIntegrityAlg>  Factory;
const IAlgFactory& EventIntegrityAlgFactory = Factory;


EventIntegrityAlg::EventIntegrityAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)
{
    declareProperty("mask"     ,  m_mask=0xffffffff); // event flag mask

}

StatusCode EventIntegrityAlg::initialize() 
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    log << MSG::INFO;
    if(log.isActive()) { 
        log.stream() <<"Applying event flag mask: " 
                     <<  std::setbase(16) <<m_mask;
    }
    log << endreq;


    return sc;
}

StatusCode EventIntegrityAlg::execute() 
{

    // Purpose and Method: Check the event flag bits and filter events as 
    //                     requested 

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

   // Retrieve the Event data for this event
 SmartDataPtr<Event::EventHeader> evtTds(eventSvc(), EventModel::EventHeader);
 if (!evtTds) {
    log << MSG::DEBUG << "No Event Header Found" << endreq;
    return sc;
 }

 SmartDataPtr<LdfEvent::EventSummaryData> summary(eventSvc(), "/Event/EventSummary");
    if( summary==0 ) {
        log << MSG::DEBUG << "No EventSummary found" << endreq;
        // if no EventSummary is available (and hence no event flag)
        // just accept the event - this might be a simulated event.. 
        // regardless we should carry on
        return sc;
     }

    unsigned int flags = summary->eventFlags();

    // If the Event Flags do exist on the TDS - check them
    // dump the event if any bit in the mask is set
    if( (m_mask!=0) && ( flags & m_mask) ) {
        // Ignoring TkrRecon Error bit
        if ( summary->badLdfStatus()) {
            setFilterPassed(false);
            log << MSG::INFO << "Event Flag contains Bad LDF Status skipping " 
                             << evtTds->event() << endreq;
         }
        if ( (summary->packetError()) || (summary->temError()) ) {
            setFilterPassed( false );
            log << MSG::INFO << "Event Flag contains Error bits - skipping " 
                             << evtTds->event() << endreq;
        }
        if (summary->badEventSeq()) {
            setFilterPassed(false);
            log << MSG::INFO << "Event Flag contains Bad Event Seq - skipping " 
                             << evtTds->event() << endreq;
        }
        if (summary->trgParityError()) {
            setFilterPassed(false);
            log << MSG::INFO << "Trigger Parity Error bit set - skipping "
                             << summary->trgParityError() << endreq;
        }
        if (summary->temBug()) {
            setFilterPassed(false);
            log << MSG::INFO << "TEM bug set - skipping "
                             << summary->temBug() << endreq;
        
        }
      
    } 

    return sc;
}

StatusCode EventIntegrityAlg::finalize() {

    setFinalized();
    StatusCode  sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    return StatusCode::SUCCESS;
}

