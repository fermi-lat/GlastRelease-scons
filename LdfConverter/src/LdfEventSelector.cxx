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

// Instantiation of a static factory class used by clients to create
// instances of this service
static const SvcFactory<LdfEventSelector> s_factory;
const ISvcFactory& LdfEventSelectorFactory = s_factory;

StatusCode LdfEventSelector::initialize()     {
    MsgStream log(msgSvc(), name());
    StatusCode sc = Service::initialize();
    if( sc.isSuccess() ) {
        setProperties();
    } else {
        log << MSG::ERROR << "Unable to initialize service " << endreq;
    }
    
    if ((m_inputList.value()).empty()) {
        log << MSG::ERROR << "No input files specified - terminating job"
            << endreq;
        return StatusCode::FAILURE;
    }

    if(m_storageType.value() != " "){
        sc = setCriteria(m_storageType);
        if(sc.isFailure()){
            log << MSG::ERROR << "Unable to parse string of input files" << endreq;
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

LdfEventSelector::LdfEventSelector( const std::string& name, ISvcLocator* svcloc ) 
: Service( name, svcloc)
, m_evtEnd(this,-1, 0) 
{

    declareProperty("StorageType", m_storageType);
    declareProperty("InputList", m_inputList);
    declareProperty( "EvtMax", m_evtMax);
    declareProperty("Instrument", m_instrument="LAT");
    declareProperty("EbfDebugLevel", m_ebfDebugLevel = 0);
    declareProperty("SweepEventSearch", m_sweepSearch = 1);
    declareProperty("GemZeroCheck", m_gemCheck=0);
    declareProperty("StartEventIndex", m_startEventIndex = 0);
    declareProperty("StartEventNumber", m_startEventNumber = 0);

    m_inputDataList = new ListName; 
    m_it = new LdfEvtIterator(this, -1, m_inputDataList->begin());
    
    //Here we get the maxEvt number from the aplication mgr property;
    //Sets the environment variable m_evtMax;
    getMaxEvent();
}


LdfEventSelector::~LdfEventSelector() {

}

bool fileExists(const std::string &fileName) {

    std::fstream checkFile;
    checkFile.open(fileName.c_str(), std::fstream::in);
    checkFile.close();
    if (checkFile.fail()) return false;
    return true;
}

/*
   Judging by the name setCriteria (part of IEvtSelector interface),
   we're not using this is the expected way.  "criteria" seems like
   it ought to have more to do with the contents of the the event,
   e.g. what the trigger was.  Instead we're using it to decide
   how to read the events in.
*/
StatusCode LdfEventSelector::setCriteria(const std::string& storageType) {
    StatusCode sc;
    MsgStream log(msgSvc(), name());

    if (storageType == "CCSDSFILE") {
      try {
      m_criteriaType = CCSDSFILE;
      m_fileName = (m_inputList.value())[0];
      facilities::Util::expandEnvVar(&m_fileName);

      if (!fileExists(m_fileName)) {
          log << MSG::ERROR << "Input file " << m_fileName 
              << " does not exist" << endreq;
          return(StatusCode::FAILURE);
      }
      
      log << MSG::DEBUG << "ctor DfiParser " << m_fileName << endreq;
      m_ebfParser = new ldfReader::DfiParser(m_fileName);
      log << MSG::DEBUG << "return from ctor" << endreq;
      } catch(...) {
        log << MSG::ERROR << "failed to setup DfiParser" << endreq;
        return(StatusCode::FAILURE);
      }
    } else if ( (storageType == "LDFFILE") || (storageType == "EBFFILE") ) {
      try {
          m_criteriaType = LDFFILE;

          m_fileName = (m_inputList.value())[0];

          // More files could be specified in the list.  Don't know what
          // Heather had in mind, so for now ignore any files beyond the
          // first one.

          facilities::Util::expandEnvVar(&m_fileName);
          if (!fileExists(m_fileName)) {
              log << MSG::ERROR << "Input file " << m_fileName 
                  << " does not exist" << endreq;
              return(StatusCode::FAILURE);
          }

          // 2nd arg. is true if FITS file, false if raw ebf
          m_ebfParser = new ldfReader::LdfParser(m_fileName, false, m_instrument);
          m_ebfParser->setDebug((m_ebfDebugLevel != 0) );
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

          m_fileName = (m_inputList.value())[0];

          // More files could be specified in the list.  Don't know what
          // Heather had in mind, so for now ignore any files beyond the
          // first one.

          facilities::Util::expandEnvVar(&m_fileName);
          if (!fileExists(m_fileName)) {
              log << MSG::ERROR << "Input file " << m_fileName 
                  << " does not exist" << endreq;
              return(StatusCode::FAILURE);
          }

          m_ebfParser = new ldfReader::LdfParser(m_fileName, true, m_instrument);
      } catch(LdfException &e) {
         log << MSG::ERROR << "LdfException: " << e.what() << endreq;
         return(StatusCode::FAILURE);
      } catch(...){
        log << MSG::ERROR << "failed to setup LdfParser" << endreq;
        return(StatusCode::FAILURE);
      }
      m_ebfParser->setDebug((m_ebfDebugLevel != 0));

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


StatusCode LdfEventSelector::setCriteria( const SelectionCriteria& ) {
    return StatusCode::SUCCESS;
}



IEvtSelector::Iterator* LdfEventSelector::begin() const {
    // Purpose and Method:  Called by ApplicationMgr::initialize( )
    MsgStream log(msgSvc(), name());
    //    StatusCode sc;       [unused for now]

    if ((m_criteriaType == CCSDSFILE) || 
        (m_criteriaType == LDFFILE) || 
        (m_criteriaType == LDFFITS) ) {
        
        log << MSG::DEBUG << " Input data set is " << m_fileName << endreq;
        m_it->m_recId = 0;                           
        
        (*m_it)++;  // increment to the first event
        
        return m_it;

    } else if(m_criteriaType == NONE) {
        
        log << MSG::DEBUG << "Using Simple counter GlastEventSelector" << endreq;
        log << MSG::INFO << "No input file selected for GlastEventSelector" << endreq;
        (*m_it)++;
        return m_it;
    } 

    log << MSG::DEBUG << "Using Simple counter LdfEventSelector" << endreq;
    log << MSG::INFO << "No input file selected for LdfEventSelector" << endreq;
    (*m_it)++;
    return m_it;
    
}


IEvtSelector::Iterator& LdfEventSelector::next(IEvtSelector::Iterator& it) 
  const {
    static bool lastEventFlag = false;
    MsgStream log(msgSvc(), name());

    log << MSG::DEBUG << "next" << endreq;
    // static counter for use when we want to skip to event N
    static unsigned int counter = 0;
    
    LdfEvtIterator* irfIt = dynamic_cast<LdfEvtIterator*>(&it);
    unsigned marker;
    try {

    if ( (lastEventFlag) && (m_criteriaType == CCSDSFILE) ) {
        // We've already found the last event, terminate loop
        *(irfIt) = m_evtEnd;
        return *irfIt;
    } 


    if  ((m_criteriaType == CCSDSFILE) ||
         (m_criteriaType == LDFFILE) ||
         (m_criteriaType == LDFFITS) )
    {
        
        irfIt->m_evtCount++;
        counter++;  //increment skip counter too

        if ((m_startEventIndex > 0) && (counter >= m_startEventIndex))  
            log << MSG::DEBUG << "Processing Event Absolute Index " << counter << endreq;
        else if ((m_startEventNumber == 0) && (m_startEventIndex == 0) )
            log << MSG::DEBUG << "Processing Event " <<  irfIt->m_evtCount << endreq;
        else if (m_startEventIndex > 0)
            log << MSG::DEBUG << "Starting Processing on Index " << m_startEventIndex << endreq;

        static bool findFirstMarkerFive = false;
        static bool findSecondMarkerFive = false;
        // Allows JO to skip the search for the sweep events
        // no skipping for LSF (CCSDSFILEs)
        if ((m_sweepSearch == 0) && (m_criteriaType != CCSDSFILE) ){
            log << MSG::WARNING << "Skipping check for first Sweep Event - ARE YOU SURE YOU WANT TO DO THIS???" << endreq;
            findFirstMarkerFive = true;
        }
        if (m_criteriaType == CCSDSFILE)
            findFirstMarkerFive = true;

        if ((m_startEventIndex > 0) && (!findFirstMarkerFive)) {
            log << MSG::WARNING << "Since StartEventIndex is non-zero, skipping search for first sweep event" << endreq;
            findFirstMarkerFive = true;
        }
        
        bool DONE=false;
        static bool foundEventNumber = false;  
        // Once we find the eventNumber, we carry on from there

        while ((!DONE) || (!findFirstMarkerFive)) {
          int status = m_ebfParser->loadData();
          if (status < 0) {
            log << MSG::INFO << "Failed to get Event" << endreq;
            log << MSG::INFO << "This job will terminate after reading" 
                << endreq;
            *(irfIt) = m_evtEnd;
            break;
          } else {

            // If a StartEventNumber JO parameter is specified, check the
            // EventSequence in the OSW contribution to see if we've reached
            // The event in question.
            if (m_startEventNumber > 0) {
                unsigned int eventSeq = ldfReader::LatData::instance()->getOsw().evtSequence();
                if (eventSeq >= m_startEventNumber) {
                    log << MSG::DEBUG << "Processing Event Number: " << eventSeq
                        << " and was searching for " << m_startEventNumber << endreq;
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
                *(irfIt) = m_evtEnd;
                return *irfIt;
            } else if (!findFirstMarkerFive) {
                if (marker == 5) { 
                    if (skippedEvents > 0)
                        log << MSG::WARNING << "Skipped " << skippedEvents 
                          << " Events before finding first sweep event" << endreq;
                    findFirstMarkerFive = true; 
                } else {
                    if (!skippedEvents)
                        log << MSG::WARNING << "First Event is not sweep event" 
                        << " marker != 5, skipping events until first "
                        << "sweep event is found" << endreq;
                    ++skippedEvents;
                }
            }

            //if ((marker != 0) || (!findFirstMarkerFive) || (counter < m_startEventIndex) || (!foundEventNumber) ) {

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
                  *(irfIt) = m_evtEnd;
                  return *irfIt;
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

        if(irfIt->m_evtCount > m_evtMax) {
            *(irfIt) = m_evtEnd;
            log << MSG::INFO << "Stopping loop at user maxEvent Request" 
                << endreq;
            if ( (m_criteriaType != CCSDSFILE) && (marker != 5) )
                log << MSG::WARNING << "Last Event was not a sweep event with"
                    << " marker == 5" << endreq;
        }
        // Move file pointer for the next event
        int ret = m_ebfParser->nextEvent();
        if ( (ret !=0 ) && (m_criteriaType == CCSDSFILE) ) { 
          log << MSG::INFO << "Last event found: source exhausted" << endreq;
          lastEventFlag = true;
        } else if (ret != 0) {
          log << MSG::INFO << "Input event source exhausted" << endreq;
          if (marker != 5) 
              log << MSG::WARNING << "Last Event was not a sweep event with"
                  << " marker == 5" << endreq;
          *(irfIt) = m_evtEnd;
        }
        return *irfIt;
    } 
    else if( m_criteriaType == NONE)
    {
        LdfEvtIterator* simpleIt = dynamic_cast<LdfEvtIterator*>(&it);
        //TODO: What the hell this really does?
        (simpleIt->m_recId)++;
        log << MSG::DEBUG << "Reading Event " <<  simpleIt->m_evtCount << endreq;
        simpleIt->m_evtCount++;
        log << MSG::DEBUG << "Reading Event " <<  simpleIt->m_evtCount << endreq;
        
        //If we go over the count set equal to the end
        if(simpleIt->m_evtCount > m_evtMax) {
            *(simpleIt) = m_evtEnd;
        }
        return *simpleIt;
    } else {
        //We are going to blow up if this situation happens
        //TODO: Error needs to be handled better
        log << MSG::FATAL << "Cannot iterate a dead iterator" << endreq;
        return it;
    }
   } catch(LdfException &e) {

      log << MSG::ERROR << "LdfException caught " << e.what() << endreq;
      *(irfIt) = m_evtEnd;
      return *irfIt;
   } catch(...) {

      log << MSG::ERROR << "Exception caught " << endreq;
      *(irfIt) = m_evtEnd;
      return *irfIt;
   }
}

IEvtSelector::Iterator* LdfEventSelector::end() const {
    IEvtSelector::Iterator* it = (IEvtSelector::Iterator*)(&m_evtEnd);
    return it;
    
}

IEvtSelector::Iterator& LdfEventSelector::previous(IEvtSelector::Iterator& it) const {
    MsgStream log(msgSvc(), name());
    log << MSG::FATAL << " LdfEventSelector Iterator, operator -- not supported " << endreq;
    
    return it;
}



IOpaqueAddress* LdfEventSelector::reference(const IEvtSelector::Iterator& it) const 
{
    // Purpose and Method:  Called from ApplicationMgr_i::nextEvent()       
    //   Create root address and assign address to data service
    //   IOpaqueAddress* addr = **m_evtIterator;
    // TO DO: understand what I need to do to in the simple case

    MsgStream log(msgSvc(), name());
    
    // convert to our iterator
    const LdfEvtIterator* irfIt = dynamic_cast<const LdfEvtIterator*>(&it);
    
    // get the file name and record id
    int recId = irfIt->m_recId;
    recId++;
    
    IOpaqueAddress* pAddr = 0;
    std::string str("");
    unsigned long temp = 0;
    if ( m_addrCreator->createAddress(TEST_StorageType, m_rootCLID, 
                                      &str, &temp, pAddr).isSuccess() ) 
    {
        return pAddr;
    }
    
    return 0;
}


StatusCode LdfEventSelector::queryInterface(const IID& riid, void** ppvInterface)  {
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

