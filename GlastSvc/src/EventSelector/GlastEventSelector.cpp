// File and Version Information:
// $Header$
// 
// Description:


#include "GlastEventSelector.h"

#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/Property.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "GaudiKernel/ObjectVector.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static const SvcFactory<GlastEventSelector> s_factory;
const ISvcFactory& GlastEventSelectorFactory = s_factory;

StatusCode GlastEventSelector::initialize()     {
    MsgStream log(msgSvc(), name());
    StatusCode sc = Service::initialize();
    if( sc.isSuccess() ) {
        setProperties();
    } else {
        log << MSG::ERROR << "Unable to initialize service " << endreq;
    }
    
    if(m_jobInput != " "){
        sc = setCriteria(m_jobInput);
        if(sc.isFailure()){
            log << MSG::ERROR << "Unable to parse string of input files" << endreq;
            return sc;
        }
    }
    // now try to find the GlastDevSvc service
    IGlastDetSvc* m_detSvc = 0;
    
    // Retrive conversion service handling event iteration
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

GlastEventSelector::GlastEventSelector( const std::string& name, ISvcLocator* svcloc ) 
: Service( name, svcloc)
, m_evtEnd(this,-1, 0) 
{
    //declareProperty( "JobInput", m_jobInput=" ");
    declareProperty( "Input", m_jobInput=" ");
    declareProperty( "EvtMax", m_evtMax);
    m_inputDataList = new ListName; 
    m_it = new GlastEvtIterator(this, -1, m_inputDataList->begin());
    // GlastEvtIterator m_evtEnd(this, -1, m_inputDataList->end());
    
    //Here we get the maxEvt number from the aplication mgr property;
    //Sets the environment variable m_evtMax;
    getMaxEvent();
}


GlastEventSelector::~GlastEventSelector() {
    delete m_inputDataList;
}

StatusCode GlastEventSelector::setCriteria( const std::string& criteria ) {
    // Purpose and Method:  Create list of input files/Job Id's

    StatusCode sc;
    MsgStream log(msgSvc(), name());
    m_criteria = criteria;
    
    // It's now time to parse the criteria for the event selection
    // The format for the criteria is:
    //        IRFFILE  filename1, filename2 ...
    std::string rest = criteria;
    std::string type;
    int ipos  = rest.find_first_not_of(' ');
    rest      = rest.substr(ipos, -1);            // remove blanks before
    int lpos  = rest.find_first_of(' ');          // locate next blank
    type      = rest.substr(0,lpos );
    
    
    if( type == "IRFFILE" ) {
        log << MSG::ERROR << "IRFFILE is not supported" << endreq;
        return StatusCode::FAILURE;
        rest      = rest.substr(lpos, -1);
    }else if( type == "NONE"){
        m_criteriaType = NONE;
        log << MSG::INFO << "GlastEventSelector input has been set to NONE" << endreq;
    }
    else {
        log << MSG::ERROR << "Invalid Event Selection Criteria: " << criteria << endreq;
        return StatusCode::FAILURE;
    }
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Invalid Event Selection Criteria: " << criteria << endreq;
        return sc;
    }
    
    IService *isvc = 0;
    sc = service("GlastDetSvc", isvc, true);
    if (sc.isSuccess() ) {
        sc = isvc->queryInterface(IID_IGlastDetSvc, (void**)&m_detSvc);
    }
    
    if(sc.isFailure()){
        log << MSG::ERROR << "Unable start Glast detector service" << endreq;
        return sc;
    }
    
    return sc;
}

StatusCode GlastEventSelector::parseStringInList( const std::string& namelist, 
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

StatusCode GlastEventSelector::setCriteria( const SelectionCriteria& criteria ) {
    return StatusCode::SUCCESS;
}


StatusCode GlastEventSelector::getFileName(ListName::const_iterator* inputIt, 
                                           std::string& fName) const {
    // Purpose and Method:  Find out the name of the file from list of files 
    //  or Jobs

    MsgStream log(msgSvc(), name());
    
    if( m_criteriaType == IRFFILE){                 // If CRITERIA = FILE Get File name
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

IEvtSelector::Iterator* GlastEventSelector::begin() const {
    // Purpose and Method:  Called by ApplicationMgr::initialize( )
    MsgStream log(msgSvc(), name());
    StatusCode sc;
    
#if 0
    if(m_criteriaType == IRFFILE){
        
        if( (m_it->m_inputDataIt) ==  m_inputDataList->end()) { 
            return m_it;                               // If no data available  return begin = end
        }
        log << MSG::DEBUG << " First input data set is " << *(m_it->m_inputDataIt) << endreq;
        m_it->m_recId = 0;                            // If data file is found Initialize record iterator
        
        std::string fName;
        sc = getFileName(&(m_it->m_inputDataIt), fName);
        
        // assumes GlastDetectors have been init. Open for read
        //it->m_inputFile = new ResponseFile (fName.c_str(), true);
        
        sc = m_detSvc->openIRF(fName);
        if(sc.isFailure()){
            *(m_it) = m_evtEnd; 
            log << MSG::INFO << "failed to open file " << fName << endreq;
            return m_it;
        }
        
        
        (*m_it)++;  // increment to the first event
        
        return m_it;
    } else if(m_criteriaType == NONE)
    {
        
        log << MSG::DEBUG << "Using Simple counter GlastEventSelector" << endreq;
        log << MSG::INFO << "No input file selected for GlastEventSelector" << endreq;
        (*m_it)++;
        return m_it;
    } 
    return 0;
#else
    log << MSG::DEBUG << "Using Simple counter GlastEventSelector" << endreq;
    log << MSG::INFO << "No input file selected for GlastEventSelector" << endreq;
    (*m_it)++;
    return m_it;
    
#endif
}


IEvtSelector::Iterator& GlastEventSelector::next(IEvtSelector::Iterator& it) const {
    MsgStream log(msgSvc(), name());
    
    if(m_criteriaType == IRFFILE)
    {
        GlastEvtIterator* irfIt = dynamic_cast<GlastEvtIterator*>(&it);
#if 0
        
        log << MSG::DEBUG << "Processing Event " <<  irfIt->m_recId << endreq;
        
        
        irfIt->m_evtCount++;
        
        // causes data from the current event in the IRF file to be loaded into the GlastDetector
        // objects
        StatusCode sc = m_detSvc->readIRF();// Ian Note: This is where the events are read in.
        
        if (sc.isFailure()){
            log << MSG::ERROR << " Failed to get Event " << irfIt->m_recId << endreq;
            log << MSG::ERROR << "The job will STOP normally after reading : " 
                << --(irfIt->m_recId) << " Events" << endreq;
            *(irfIt) = m_evtEnd;
        }
        
        if(irfIt->m_evtCount > m_evtMax) {
            *(irfIt) = m_evtEnd;
            log << MSG::INFO << "Stoping loop at user maxEvent Request" << endreq;
        }
#endif
        return *irfIt;
    } 
    else if( m_criteriaType == NONE)
    {
        GlastEvtIterator* simpleIt = dynamic_cast<GlastEvtIterator*>(&it);
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

IEvtSelector::Iterator* GlastEventSelector::end() const {
    IEvtSelector::Iterator* it = (IEvtSelector::Iterator*)(&m_evtEnd);
    return it;
    
}

IEvtSelector::Iterator& GlastEventSelector::previous(IEvtSelector::Iterator& it) const {
    MsgStream log(msgSvc(), name());
    log << MSG::FATAL << " GlastEventSelector Iterator, operator -- not supported " << endreq;
    
    return it;
}



IOpaqueAddress* GlastEventSelector::reference(const IEvtSelector::Iterator& it) const 
{
    // Purpose and Method:  Called from ApplicationMgr_i::nextEvent()       
    //   Create root address and assign address to data service
    //   IOpaqueAddress* addr = **m_evtIterator;
    // TO DO: understand what I need to do to in the simple case

    MsgStream log(msgSvc(), name());
    
    // convert to our iterator
    const GlastEvtIterator* irfIt = dynamic_cast<const GlastEvtIterator*>(&it);
    
    // get the file name and record id
    const std::string& fName = *(irfIt->m_inputDataIt);
    int recId = irfIt->m_recId;
    recId++;
    
    ListName::const_iterator* inputDataIt = (ListName::const_iterator*)(&irfIt->m_inputDataIt);
    
    
    IOpaqueAddress* pAddr = 0;
    std::string str("");
    unsigned long temp = 0;
    if ( m_addrCreator->createAddress(SICB_StorageType, m_rootCLID, &str, &temp, pAddr).isSuccess() ) 
    {
        return pAddr;
    }
    
    return 0;
}


StatusCode GlastEventSelector::queryInterface(const IID& riid, void** ppvInterface)  {
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

StatusCode GlastEventSelector::getMaxEvent()
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

