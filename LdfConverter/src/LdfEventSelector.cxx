// File and Version Information:
// $Header$
// 
// Description:


#include "LdfEventSelector.h"

#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/ObjectVector.h"

#include "ldfReader/data/LatData.h"
#include "ldfReader/LdfException.h"
#include <exception>
#include <fstream>
#include "facilities/Util.h"
#include "ldfReader/LdfParser.h"
#include "ldfReader/DfiParser.h"
#include "ldfReader/SocketParser.h"
#include "LdfSelectorContext.h"

#include <sys/stat.h>

// Instantiation of a static factory class used by clients to create
// instances of this service
static const SvcFactory<LdfEventSelector> s_factory;
const ISvcFactory& LdfEventSelectorFactory = s_factory;


///NEW

/// Not currently supported for LDF/LSF reading
StatusCode LdfEventSelector::resetCriteria(const std::string& /* criteria */,
                                           Context& /* context  */)  const
{
    return StatusCode::FAILURE;
}


/// Create a new event loop context
StatusCode LdfEventSelector::createContext(Context*& refpCtxt) const
{
    // Max event is zero. Return begin = end
    refpCtxt = 0;
    if ( m_firstEvent < 0 ) {
        MsgStream log(messageService(), name());
        log << MSG::ERROR  << "First Event = " << m_firstEvent << " not valid" 
            << endreq;
        log << MSG::ERROR  << "It should be > 0 " << endreq;
        return StatusCode::FAILURE;    // if failure => iterators = end();
    }
    LdfSelectorContext* ctxt = new LdfSelectorContext(this);
    //HMKctxt->set(0, -1, 0, 0);
    ctxt->set(0,0,ctxt,0);
    refpCtxt = ctxt;
    long nskip = m_firstEvent;
    while( --nskip > 0 )    {
        StatusCode sc = next(*refpCtxt);
        if ( sc.isFailure() ) {
            MsgStream log(messageService(), name());
            log << MSG::ERROR 
                << " createContext() failed to start with event number "
                << m_firstEvent << endreq;
            releaseContext(refpCtxt);
            refpCtxt = 0;
            return StatusCode::FAILURE;
        }
    }
    return StatusCode::SUCCESS;
}

/// Get next iteration item from the event loop context
StatusCode LdfEventSelector::next(Context& refCtxt) const  {
    return next(refCtxt, 1);
}


/// Get previous iteration item from the event loop context
StatusCode LdfEventSelector::previous(Context& refCtxt) const  {
    return previous(refCtxt, 1);
}

/// Get previous iteration item from the event loop context, but skip jump elements
StatusCode LdfEventSelector::previous(Context& it,int jump) const  {
    // HMK NOT SUPPORTED FOR LDF YET
    /*
    LdfSelectorContext *pIt  = dynamic_cast<LdfSelectorContext*>(&it);
    if ( pIt && jump > 0 )    {
    StatusCode sc = StatusCode::SUCCESS;
    for ( int i = 0; i < jump && sc.isSuccess(); ++i ) {
    const LdfEventSelectorDataStream* s = m_streams[pIt->ID()];
    Context* it = pIt->context();
    IEvtSelector* sel = s->selector();
    if ( it && sel )    { // First explot the current stream
    // This stream is empty: advance to the next stream
    sc = sel->previous(*it);  // This stream is empty: advance to the next stream
    if ( !sc.isSuccess() )   {
    sc = lastOfPreviousStream(true, *pIt);
    }
    else  {
    pIt->increaseCounters(false);
    pIt->set(it, 0);
    }
    printEvtInfo(pIt);
    if ( !sc.isSuccess() ) {
    return sc;
    }
    }
    pIt->increaseCounters(false);
    }
    return sc;
    }
    printEvtInfo(pIt);
    */
    return StatusCode::FAILURE;
}

/// Access last item in the iteration - not yet supported
StatusCode LdfEventSelector::last(Context& refCtxt) const {
    LdfSelectorContext *pIt  = dynamic_cast<LdfSelectorContext*>(&refCtxt);
    if ( pIt )    {
    }
    return StatusCode::FAILURE;
}

/// Rewind the dataset
StatusCode LdfEventSelector::rewind(Context& refCtxt) const  {
    // NOT SUPPORTED IN LDF YET HMK
    /*
    LdfSelectorContext *ctxt  = dynamic_cast<LdfSelectorContext*>(&refCtxt);
    if ( ctxt )    {
    ctxt->set(0, -1, 0, 0);
    firstOfNextStream(true, *ctxt);
    long nskip = m_firstEvent;
    while( --nskip > 0 )    {
    StatusCode sc = next(*ctxt);
    if ( sc.isFailure() ) {
    MsgStream log(messageService(), name());
    log << MSG::ERROR << "rewind() failed to start with event number "
    << m_firstEvent << endreq;
    return StatusCode::FAILURE;
    }
    }
    return StatusCode::SUCCESS;
    } 
    */
    return StatusCode::FAILURE;
}

/// Create new Opaque address corresponding to the current record
StatusCode LdfEventSelector::createAddress(const Context&   refCtxt,
                                           IOpaqueAddress*& refpAddr) const
{
    const LdfSelectorContext *cpIt  = dynamic_cast<
        const LdfSelectorContext*>(&refCtxt);
    LdfSelectorContext *pIt  = const_cast<LdfSelectorContext*>(cpIt);
    refpAddr = 0;
    IOpaqueAddress* pAddress = 0;
    std::string str("");
    unsigned long temp = 0;
    if ( m_addrCreator->createAddress(TEST_StorageType, m_rootCLID, &str, &temp, pAddress).isSuccess() ) 
    {
        refpAddr = pAddress;
        return StatusCode::SUCCESS;
    } else    {
        MsgStream log(msgSvc(), name());
        log << MSG::ERROR << "Failed to createAddress " 
            << " CLID:"   << m_rootCLID
            << endreq;
        return StatusCode::FAILURE;
    }

    return StatusCode::FAILURE;

}

// Release existing event iteration context
StatusCode LdfEventSelector::releaseContext(Context*& refCtxt) const  {
    const LdfSelectorContext *cpIt = dynamic_cast<
        const LdfSelectorContext*>(refCtxt);
    LdfSelectorContext       *pIt  = const_cast<LdfSelectorContext*>(cpIt);
    if ( pIt && pIt->ID() >= 0 )    {
       refCtxt = 0;
       delete pIt;
       return StatusCode::SUCCESS;
    }
    return StatusCode::FAILURE;
}


// Progress report
void LdfEventSelector::printEvtInfo(const LdfSelectorContext* iter) const {
    if ( 0 != iter )  {
        long count = iter->numEvent();
        // Print an message every m_evtPrintFrequency events
        if ( 0 == iter->context() )   {
            MsgStream log(messageService(), name());
            log << MSG::INFO << "End of event input reached." << endreq;
        }
        else if( iter->numStreamEvent() == -1 ) {
            // Intial value for this stream
        }
        else if( m_evtPrintFrequency != -1 && (count % m_evtPrintFrequency == 0))
        {
            MsgStream log(messageService(), name());
            log << MSG::ALWAYS << "Reading Event record " << count+1
                << ". Record number within stream " << iter->ID()+1
                << ": " << iter->numStreamEvent()+1 << endreq;
        }
    }
    else  {
        MsgStream log(messageService(), name());
        log << MSG::INFO << "End of event input reached." << endreq;
    }
} 


/// OLD




StatusCode LdfEventSelector::initialize()     {
    MsgStream log(msgSvc(), name());
    StatusCode sc = Service::initialize();
    if( sc.isSuccess() ) {
        setProperties();
    } else {
        log << MSG::ERROR << "Unable to initialize service " << endreq;
    }

    if ( (m_fileName.empty()) && (m_socket != 1) ) {
        log << MSG::ERROR << "No input file specified - terminating job"
            << endreq;
        return StatusCode::FAILURE;
    }

    if(m_storageType.value() != " "){
        sc = setCriteria(m_storageType);
        if(sc.isFailure()){
            log << MSG::ERROR << "Unable to parse string of input files" 
                << endreq;
            return sc;
        }
    } else {
        log << MSG::ERROR << "No Storage Type specified - terminating job"
            << endreq;
        return StatusCode::FAILURE;
    }

    // Retrieve conversion service handling event iteration
    sc = serviceLocator()->service("EventCnvSvc", m_addrCreator);
    if( !sc.isSuccess() ) {
        log << MSG::ERROR << 
            "Unable to localize interface IID_IAddressCreator from service:" 
            << "EventCnvSvc" 
            << endreq;
        return sc;
    }
    // Get DataSvc
    IDataManagerSvc* eds = 0;
    sc = serviceLocator()->service("EventDataSvc", eds, true);
    if( !sc.isSuccess() ) {
        log << MSG::ERROR 
            << "Unable to localize interface IID_IDataManagerSvc "
            << "from service EventDataSvc"
            << endreq;
        return sc;
    }

    m_rootCLID = eds->rootCLID();

    return sc;
}

LdfEventSelector::LdfEventSelector( const std::string& name, 
                                   ISvcLocator* svcloc ) 
                                   : Service( name, svcloc)
{
    // From new Gaudi v18r1 EventSelector in GaudiSvc
    m_firstEvent        = 0;
    m_evtPrintFrequency = 10;
    m_evtMax            = INT_MAX;  // Obsolete parameter
    declareProperty( "FirstEvent", m_firstEvent);
    declareProperty( "EvtMax",     m_evtMax);
    declareProperty( "PrintFreq",  m_evtPrintFrequency);
    m_reconfigure = false;

    declareProperty("StorageType", m_storageType);
    declareProperty("Instrument", m_instrument="LAT");
    declareProperty("EbfDebugLevel", m_ebfDebugLevel = 0);
    declareProperty("SweepEventSearch", m_sweepSearch = 1);
    declareProperty("GemZeroCheck", m_gemCheck=0);
    declareProperty("StartEventIndex", m_startEventIndex = 0);
    declareProperty("StartEventNumber", m_startEventNumber = 0);
    declareProperty("FileName", m_fileName="");
    declareProperty("AcdRemapFile", m_acdRemap="");
    declareProperty("IgnoreSegFault", m_ignoreSegFault=0);
    declareProperty("OldStyleRunId", m_oldStyleRunId=0);

    // Options for using socket connections
    declareProperty("SocketConnection", m_socket = 0);
    // default server is 60 - chosen for the beamtest
    declareProperty("SocketServer", m_server = 60);

    //Here we get the maxEvt number from the aplication mgr property;
    //Sets the environment variable m_evtMax;
    getMaxEvent();
}


LdfEventSelector::~LdfEventSelector() {

}

bool fileExists(const std::string &fileName) {

#ifdef WIN32
    struct _stat checkFile;
    if (_stat(fileName.c_str(), &checkFile) == -1)
#else
    struct stat checkFile;
    if (stat(fileName.c_str(), &checkFile) == -1)
#endif
      return false;
    return true;

    //std::fstream checkFile;
    //checkFile.open(fileName.c_str(), std::fstream::in);
    //checkFile.close();
    //if (checkFile.fail()) return false;
    //return true;
}

/*
Judging by the name setCriteria (part of IEvtSelector interface),
we're not using this is the expected way.  "criteria" seems like
it ought to have more to do with the contents of the the event,
e.g. what the trigger was.  Instead we're using it to decide
how to read the events in.
*/
StatusCode LdfEventSelector::setCriteria(const std::string& storageType) {
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Check if using socket connections first
    if (m_socket) {
        try {
            // sockets assume LDF files
            m_criteriaType = LDFFILE;

            m_ebfParser = new ldfReader::SocketParser(m_server);
            m_ebfParser->setDebug(m_ebfDebugLevel);
            m_ebfParser->setOldStyleRunId((m_oldStyleRunId != 0));
            if (m_acdRemap != "")
                if (m_ebfParser->setAcdRemap(m_acdRemap) < 0) 
                    log << MSG::WARNING << "Failed to read ACD remap file" 
                        << endreq;
                else 
                    log << MSG::INFO << "Opened ACD remap file: " <<
                           m_acdRemap << endreq;
            if (m_ignoreSegFault)
                m_ebfParser->setIgnoreSegFault();
            return(StatusCode::SUCCESS);
        } catch(...) {
            log << MSG::ERROR << "Failed to set up socket connection" << endreq;
            return(StatusCode::FAILURE);
        }
    }

    if (storageType == "CCSDSFILE") {
        try {
            m_criteriaType = CCSDSFILE;
            facilities::Util::expandEnvVar(&m_fileName);

            if (!fileExists(m_fileName)) {
                log << MSG::ERROR << "Input file " << m_fileName 
                    << " does not exist" << endreq;
                return(StatusCode::FAILURE);
            }

            log << MSG::DEBUG << "ctor DfiParser " << m_fileName << endreq;
            m_ebfParser = new ldfReader::DfiParser(m_fileName);
            m_ebfParser->printHeader();
            m_ebfParser->setDebug(m_ebfDebugLevel);
            m_ebfParser->setOldStyleRunId((m_oldStyleRunId != 0));
            log << MSG::DEBUG << "return from ctor" << endreq;
        } catch(...) {
            log << MSG::ERROR << "failed to setup DfiParser" << endreq;
            return(StatusCode::FAILURE);
        }
    } else if ( (storageType == "LDFFILE") || (storageType == "EBFFILE") ) {
        try {
            m_criteriaType = LDFFILE;

            facilities::Util::expandEnvVar(&m_fileName);
            if (!fileExists(m_fileName)) {
                log << MSG::ERROR << "Input file " << m_fileName 
                    << " does not exist" << endreq;
                return(StatusCode::FAILURE);
            }

            log << MSG::DEBUG << "Starting LdfParser with " 
                << m_fileName << endreq;

            // 2nd arg. is true if FITS file, false if raw ebf
            m_ebfParser = new ldfReader::LdfParser(m_fileName, false, m_instrument);
            if (!m_ebfParser) 
                log << MSG::DEBUG << "Creating LdfParser failed" << endreq;
            m_ebfParser->setDebug(m_ebfDebugLevel);
            m_ebfParser->setOldStyleRunId((m_oldStyleRunId != 0));
            if (m_acdRemap != "")
                if (m_ebfParser->setAcdRemap(m_acdRemap) < 0) 
                    log << MSG::WARNING << "Failed to read ACD remap file" << endreq;
                else 
                    log << MSG::INFO << "Opened ACD remap file: " <<
                           m_acdRemap << endreq;
            if (m_ignoreSegFault)
                m_ebfParser->setIgnoreSegFault();

        } catch(LdfException &e) {
            log << MSG::ERROR << "LdfException: " << e.what() << endreq;
            return(StatusCode::FAILURE);
        } catch(...) {
            log << MSG::ERROR << "failed to setup LdfParser" << endreq;
            return(StatusCode::FAILURE);
        }

    } else if ((storageType == "LDFFITS") || (storageType == "EBFFITS")) {
        try {
            m_criteriaType = LDFFITS;

            facilities::Util::expandEnvVar(&m_fileName);
            if (!fileExists(m_fileName)) {
                log << MSG::ERROR << "Input file " << m_fileName 
                    << " does not exist" << endreq;
                return(StatusCode::FAILURE);
            }

            m_ebfParser = new ldfReader::LdfParser(m_fileName, true, m_instrument);
            m_ebfParser->setDebug(m_ebfDebugLevel);
            m_ebfParser->setOldStyleRunId((m_oldStyleRunId != 0));
            if (m_acdRemap != "")
                if (m_ebfParser->setAcdRemap(m_acdRemap) < 0) 
                    log << MSG::WARNING << "Failed to read ACD remap file" << endreq;
                else 
                    log << MSG::INFO << "Opened ACD remap file: " <<
                           m_acdRemap << endreq;
            if (m_ignoreSegFault)
                m_ebfParser->setIgnoreSegFault();

        } catch(LdfException &e) {
            log << MSG::ERROR << "LdfException: " << e.what() << endreq;
            return(StatusCode::FAILURE);
        } catch(...){
            log << MSG::ERROR << "failed to setup LdfParser" << endreq;
            return(StatusCode::FAILURE);
        }

    }else if (storageType == "NONE"){
        m_criteriaType = NONE;
        log << MSG::INFO << "LdfEventSelector input has been set to NONE" 
            << endreq;
    }
    else {
        log << MSG::ERROR << "Invalid Event Selection Criteria: " 
            << m_criteriaType << endreq;
        return StatusCode::FAILURE;
    }
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Invalid Event Selection Criteria: " 
            << m_criteriaType << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}


/// Get next iteration item from the event loop context, but skip jump elements
StatusCode LdfEventSelector::next(Context& refCtxt, int /* jump */ ) const  {
    static bool lastEventFlag = false;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << "next" << endreq;
    // static counter for use when we want to skip to event N
    static unsigned long long counter = 0;
    unsigned marker;
    try {
        LdfSelectorContext *pIt  = dynamic_cast<LdfSelectorContext*>(&refCtxt);
        if ( pIt )    {
            if ( pIt->ID() != -1 ) {


                if ( (lastEventFlag) && (m_criteriaType == CCSDSFILE) ) {
                    // We've already found the last event, terminate loop
                    return StatusCode::FAILURE;
                } 


                if  ((m_criteriaType == CCSDSFILE) ||
                    (m_criteriaType == LDFFILE) ||
                    (m_criteriaType == LDFFITS) )
                {

                    pIt->increaseCounters(false); // do this here??
                    counter++;  //increment skip counter too

                    if ((m_startEventIndex > 0) && (counter >= m_startEventIndex))  
                        log << MSG::DEBUG << "Processing Event Absolute Index " 
                        << (long long int) counter << endreq;
                    else if ((m_startEventNumber == 0) && (m_startEventIndex == 0) )
                        log << MSG::DEBUG << "Processing Event " 
                        << pIt->numEvent() << endreq;
                    else if (m_startEventIndex > 0)
                        log << MSG::DEBUG << "Starting Processing on Index " 
                        << m_startEventIndex << endreq;

                    static bool findFirstMarkerFive = false;
                    static bool findSecondMarkerFive = false;
                    // Allows JO to skip the search for the sweep events
                    // no skipping for LSF (CCSDSFILEs)
                    if ((m_sweepSearch == 0) && (m_criteriaType != CCSDSFILE) ){
                        static short sweepFlag = 0;
                        if (sweepFlag < 5)
                            log << MSG::WARNING 
                                << "Skipping check for first Sweep Event - "
                                << " ARE YOU SURE YOU WANT TO DO THIS???" 
                                << endreq;
                        findFirstMarkerFive = true;
                        ++sweepFlag;
                    }
                    if (m_criteriaType == CCSDSFILE) findFirstMarkerFive = true;

                    if ((m_startEventIndex > 0) && (!findFirstMarkerFive)) {
                        log << MSG::WARNING 
                            << "Since StartEventIndex is non-zero, "
                            << "skipping search for first sweep event" 
                            << endreq;
                        findFirstMarkerFive = true;
                    }

                    bool DONE=false;
                    static bool foundEventNumber = false;  
                    // Once we find the eventNumber, we carry on from there

                    while ((!DONE) || (!findFirstMarkerFive)) {
                        int status = m_ebfParser->loadData();
                        if (status < 0) {
                            log << MSG::INFO << "Failed to get Event" << endreq;
                            log << MSG::INFO 
                                << "This job will terminate after reading" 
                                << endreq;
                            return StatusCode::FAILURE; // HMK thnk I return failure here for context  - that will terminate loop
                            break;
                        } else {


                            // Check for ADF Header or Trailer
                            if(ldfReader::LatData::instance()->adfHdrTlr()) {
                                // Move file pointer for the next event
                                log << MSG::DEBUG << "Found ADF Hdr or Tlr" << endreq;
                                int ret = m_ebfParser->nextEvent();
                                if (ret != 0) {
                                    log << MSG::INFO << "Input event source exhausted" << endreq;
                                    if (counter != m_ebfParser->eventCount())
                                        log << MSG::WARNING << "Number of events process " 
                                            << (long long int) counter 
                                            << " does not match number of events in input file "
                                            << (long long int) m_ebfParser->eventCount() 
                                            << endreq;
                                    else
                                       log << MSG::INFO << "Processed " 
                                           << (long long int) counter 
                                           << " event from a file containing " 
                                           << (long long int) m_ebfParser->eventCount() 
                                           << " events" << endreq;
                                    return StatusCode::FAILURE;  // not sure if we're skipping the last event or not
                                }
                                continue;
                            }

                            // If a StartEventNumber JO parameter is specified, check the
                            // EventSequence if we've reached the event in question.
                            if (m_startEventNumber > 0) {
                                unsigned long long eventSeq = ldfReader::LatData::instance()->eventId();
                                if (eventSeq >= m_startEventNumber) {
                                    log << MSG::DEBUG << "Processing Event Number: " 
                                        << (long long int) eventSeq
                                        << " and was searching for " << m_startEventNumber 
                                        << endreq;
                                    foundEventNumber = true;
                                }
                            } else 
                                foundEventNumber = true;

                            // Check marker to see if this is a data event
                            unsigned int summary = ldfReader::LatData::instance()->summaryData().summary();
                            marker = EventSummary::marker(summary);
                            static unsigned int skippedEvents = 0;

                            if ( (m_criteriaType != CCSDSFILE) && (findFirstMarkerFive) 
                                && (m_sweepSearch == 1) && (marker == 5) ) {
                                    log << MSG::WARNING << "Found second marker 5 sweep event in"
                                        << " non-LSF file - terminating run" << endreq;
                                    return StatusCode::FAILURE;
                                } else if (!findFirstMarkerFive) {
                                    if (marker == 5) { 
                                        if (skippedEvents > 0)
                                            log << MSG::WARNING << "Skipped " << skippedEvents 
                                            << " Events before finding first sweep event" 
                                            << endreq;
                                        findFirstMarkerFive = true; 
                                    } else {
                                        if (!skippedEvents)
                                            log << MSG::WARNING << "First Event is not sweep event" 
                                            << " marker != 5, skipping events until first "
                                            << "sweep event is found" << endreq;
                                        ++skippedEvents;
                                    }
                                }


                                // If we haven't found first Sweep event, or haven't located the
                                // event index we're searching for or found the event number or
                                // we're reading LDF and have marker == 5
                                // THEN read in next event and continue
                                if ( (!findFirstMarkerFive) ||
                                    (counter < m_startEventIndex) ||
                                    (!foundEventNumber) ||
                                    ( (m_criteriaType != CCSDSFILE) && (marker == 5) ) ) {

                                        // Always increment skip counter, no matter what type of data 
                                        // was loaded.
                                        counter++;

                                        // Move file pointer for the next event
                                        int ret = m_ebfParser->nextEvent();
                                        if (ret != 0) {
                                            log << MSG::INFO << "Input event source exhausted" << endreq;
                                            if (counter != m_ebfParser->eventCount())
                                                log << MSG::WARNING << "Number of events process " 
                                                << (long long int) counter 
                                                << " does not match number of events in input file "
                                                << (long long int) m_ebfParser->eventCount() 
                                                << endreq;
                                            else
                                                log << MSG::INFO << "Processed " 
                                                << (long long int) counter 
                                                << " event from a file containing " 
                                                << (long long int) m_ebfParser->eventCount() 
                                                << " events" << endreq;
                                            return StatusCode::FAILURE;  // not sure if we're skipping the last event or not
                                        }
                                        // Otherwise, check to see if we found an event to process
                                    } else if ( (findFirstMarkerFive) && 
                                        (counter >= m_startEventIndex) && 
                                        (foundEventNumber) ) 
                                        DONE = true;
                        }
                    }


                    if (m_gemCheck) {
                        static short foundFive = 0;
                        if ( (foundFive < 5) &&
                            (ldfReader::LatData::instance()->getGem().conditionSummary() == 0) )  {
                                if (foundFive < 4)
                                    log << MSG::WARNING << "GEM Condition Summary Zero" << endreq;
                                else if (foundFive == 4)
                                    log << MSG::WARNING << "GEM Condition Summary Zero - further warnings will be suppressed" << endreq; 
                                ++foundFive;
                            }

                    }

                    if (pIt->numEvent() > m_evtMax) {
                        log << MSG::INFO << "Stopping loop at user maxEvent Request" 
                            << endreq;
                        if ( (m_criteriaType != CCSDSFILE) && (marker != 5) )
                            log << MSG::WARNING << "Last Event was not a sweep event with"
                            << " marker == 5" << endreq;
                        return StatusCode::FAILURE;
                    }
                    // Move file pointer for the next event
                    int ret = m_ebfParser->nextEvent();
                    if ( (ret !=0 ) && (m_criteriaType == CCSDSFILE) ) { 
                        log << MSG::INFO << "Last event found: source exhausted" << endreq;
                        lastEventFlag = true;
                        if (counter != m_ebfParser->eventCount())
                            log << MSG::WARNING << "Number of events process " 
                            << (long long int) counter
                            << " does not match number of events in input file "
                            << (long long int) m_ebfParser->eventCount() << endreq;
                        else
                            log << MSG::INFO << "Processed " << (long long int) counter 
                            << " event from a file containing " 
                            << (long long int) m_ebfParser->eventCount() 
                            << " events" << endreq;
                    } else if (ret != 0) {
                        log << MSG::INFO << "Input event source exhausted" << endreq;
                        if (marker != 5) 
                            log << MSG::WARNING << "Last Event was not a sweep event with"
                            << " marker == 5" << endreq;
                        if (counter != m_ebfParser->eventCount())
                            log << MSG::WARNING << "Number of events process " 
                            << (long long int) counter
                            << " does not match number of events in input file "
                            << (long long int) m_ebfParser->eventCount() << endreq;
                        else 
                            log << MSG::INFO << "Processed " << (long long int) counter 
                            << " event from a file containing " 
                            << (long long int) m_ebfParser->eventCount() 
                            << " events" << endreq;
                        return StatusCode::FAILURE;
                        // HMK*(irfIt) = m_evtEnd;
                        //maybe need a terminate flag
                    }
                    return StatusCode::SUCCESS;
                    //HMKreturn *irfIt;
                } else {
                    //We are going to blow up if this situation happens
                    //TODO: Error needs to be handled better
                    log << MSG::FATAL << "Cannot iterate a dead iterator" << endreq;
                    return StatusCode::FAILURE;
                }
            } 
        }  //  end of context ifs
    } catch(LdfException &e) {

        log << MSG::ERROR << "LdfException caught " << e.what() << endreq;
        return StatusCode::FAILURE;
    } catch(...) {

        log << MSG::ERROR << "Exception caught " << endreq;
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}


StatusCode LdfEventSelector::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
    if ( riid == IID_IEvtSelector )  {
        *ppvInterface = (IEvtSelector*)this;
    }
    else if ( riid == IID_IProperty )  {
        *ppvInterface = (IProperty*)this;
    }
    else   {
        return Service::queryInterface( riid, ppvInterface );
    }
    addRef();
    return SUCCESS;
}

StatusCode LdfEventSelector::getMaxEvent()
{
    IProperty* appPropMgr=0;
    StatusCode status = 
        serviceLocator()->getService("ApplicationMgr", IID_IProperty,
        reinterpret_cast<IInterface*&>( appPropMgr ));
    if( status.isFailure() ) return status;

    IntegerProperty evtMax("EvtMax",0);
    status = appPropMgr->getProperty( &evtMax );
    if (status.isFailure()) return status;

    int max_event = evtMax.value();

    status = appPropMgr->setProperty( evtMax );
    appPropMgr->release();
    m_evtMax = max_event;
    return status;
}

