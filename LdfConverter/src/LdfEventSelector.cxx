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

#include "facilities/Util.h"

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
    
    if(m_storageType.value() != " "){
        sc = setCriteria(m_storageType);
        if(sc.isFailure()){
            log << MSG::ERROR << "Unable to parse string of input files" << endreq;
            return sc;
        }
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

    m_inputDataList = new ListName; 
    m_it = new LdfEvtIterator(this, -1, m_inputDataList->begin());
    
    //Here we get the maxEvt number from the aplication mgr property;
    //Sets the environment variable m_evtMax;
    getMaxEvent();
}


LdfEventSelector::~LdfEventSelector() {

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

    if (storageType == "EBFFILE") {
      m_criteriaType = EBFFILE;

      m_fileName = (m_inputList.value())[0];

      // More files could be specified in the list.  Don't know what
      // Heather had in mind, so for now ignore any files beyond the
      // first one.

      facilities::Util::expandEnvVar(&m_fileName);

      // 2nd arg. is true if FITS file, false if raw ebf
      try {
          m_ebfParser = new ldfReader::LdfParser(m_fileName, false, m_instrument);
      } catch(...) {
        log << MSG::ERROR << "failed to setup EbfParser" << endreq;
        return(StatusCode::FAILURE);
      }
      m_ebfParser->setDebug((m_ebfDebugLevel != 0) );

    } else if (storageType == "EBFFITS") {
      m_criteriaType = EBFFITS;

      m_fileName = (m_inputList.value())[0];

      // More files could be specified in the list.  Don't know what
      // Heather had in mind, so for now ignore any files beyond the
      // first one.

      facilities::Util::expandEnvVar(&m_fileName);

      try {
          m_ebfParser = new ldfReader::LdfParser(m_fileName, true, m_instrument);
      } catch(...){
        log << MSG::ERROR << "failed to setup EbfParser" << endreq;
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
        return sc;
    }
        
    return sc;
}

/*
// Doesn't look like this is needed
StatusCode LdfEventSelector::parseStringInList( const std::string& namelist, 
                                                 ListName* inputDataList ) {
    // Purpose and Method: Parse criteria string: Fill in the list of input 
    //  files or Job Id's.  Also parse out an environment variable and 
    //  substitute for it.
    if(m_criteriaType != NONE)
    {
        std::string rest = namelist;
        std::string substitute;
        while(true) {
            int ipos = rest.find_first_not_of(" ,");
            if (ipos == -1 ) break;
            rest = rest.substr(ipos, -1);            // remove blanks before
            int lpos  = rest.find_first_of(" ,");    // locate next blank
            if (lpos == -1 ) {
                rest = rest.substr(0,lpos );
                
                //now pull out and substitute for environment vaiables
                int envStart = rest.find_first_of("$(");
                int envEnd = rest.find_first_of(")");
                
                // add 2 characters to get rid of $(
                int afterBracket = envStart + 2;
                
                if(!((envStart==-1)||(envEnd==-1)))
                {
                    std::string envVariable = rest.substr(afterBracket,(envEnd-afterBracket));
                    const char * instruPath = ::getenv(envVariable.data());
                    substitute = rest.replace(envStart,(envEnd+1), instruPath);
                    inputDataList->push_back(substitute);       // insert last item in list and
                    break;
                }
                inputDataList->push_back(rest);    
                break;                                 
            }
            inputDataList->push_back( rest.substr(0,lpos ));   // insert in list
            rest = rest.substr(lpos, -1);                      // get the rest
        }
        
        return StatusCode::SUCCESS;
    } else {
        return StatusCode::FAILURE;
    }
}
*/

StatusCode LdfEventSelector::setCriteria( const SelectionCriteria& ) {
    return StatusCode::SUCCESS;
}


/*
StatusCode LdfEventSelector::getFileName(ListName::const_iterator* inputIt, 
                                           std::string& fName) const {
    // Purpose and Method:  Find out the name of the file from list of files 
    //  or Jobs

    MsgStream log(msgSvc(), name());
    
    if( m_criteriaType == EBFFILE){                 // If CRITERIA = FILE Get File name
        fName = **inputIt;
    } else if(m_criteriaType == NONE)  {
        log << MSG::ERROR << "Selection Criteria set to NONE can't get file Name" << endreq;
        return StatusCode::FAILURE;
    } else {
        log << MSG::ERROR << "Wrong Selection Criteria, either NONE or IRFFILE" << endreq;
        return StatusCode::FAILURE;
    }
    
    return StatusCode::SUCCESS;
}
*/

IEvtSelector::Iterator* LdfEventSelector::begin() const {
    // Purpose and Method:  Called by ApplicationMgr::initialize( )
    MsgStream log(msgSvc(), name());
    //    StatusCode sc;       [unused for now]

    if ((m_criteriaType == EBFFILE) || (m_criteriaType == EBFFITS) ) {
        
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
    MsgStream log(msgSvc(), name());
    
    if  ((m_criteriaType == EBFFILE) ||
         (m_criteriaType == EBFFITS) )
    {
        LdfEvtIterator* irfIt = dynamic_cast<LdfEvtIterator*>(&it);
        
        //log << MSG::DEBUG << "Processing Event " <<  irfIt->m_recId << endreq;
        log << MSG::DEBUG << "Processing Event " <<  irfIt->m_evtCount << endreq;
        
        irfIt->m_evtCount++;
        
        int status = m_ebfParser->loadData();
        if (status < 0) {
            log << MSG::ERROR << "Failed to get Event" << endreq;
            log << MSG::ERROR << "This job will terminate after reading" 
                << endreq;
            *(irfIt) = m_evtEnd;
        }
        
        if(irfIt->m_evtCount > m_evtMax) {
            *(irfIt) = m_evtEnd;
            log << MSG::INFO << "Stopping loop at user maxEvent Request" 
                << endreq;
        }
        // Move file pointer for the next event
        int ret = m_ebfParser->nextEvent();
        if (ret != 0) {
          log << MSG::INFO << "Input event source exhausted" << endreq;
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

