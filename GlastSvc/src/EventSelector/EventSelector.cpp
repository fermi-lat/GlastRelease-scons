// $Header$
//====================================================================
//  EventSelector.cpp
//--------------------------------------------------------------------
//
//  Package    : EventSelector  (The Glast IRF Event Selector Package)
//
//
//====================================================================

// Include files
#include "EventSelector.h"

#include "Gaudi/MessageSvc/MsgStream.h"

#include "Gaudi/Kernel/StatusCode.h"
#include "Gaudi/Kernel/SvcFactory.h"
#include "Gaudi/Kernel/AddrFactory.h"

#include "Gaudi/Interfaces/IDataProviderSvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "instrument/DetectorConverter.h"
#include "instrument/Scintillator.h"

#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/MonteCarlo/MCACDHit.h"
#include "GlastEvent/MonteCarlo/MCCalorimeterHit.h" //TODO: remove this declaration

extern const CLID& CLID_Event;

// Glast address factory
extern unsigned char GLAST_StorageType;
extern const IAddrFactory& DummyAddressFactory;

// Instantiation of a static factory class used by clients to create
// instances of this service
static const SvcFactory<EventSelector> s_factory;
const ISvcFactory& GlastEventSelectorFactory = s_factory;

StatusCode EventSelector::initialize()     {
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
    
	//StatusCode sc;
	sc = serviceLocator()->getService( "EventDataSvc", IID_IDataProviderSvc, (IInterface*&)m_eventDataSvc );
	if( sc.isFailure() ) {
		log << MSG::ERROR << "Unable to locate service " << "EventDataSvc" << endreq;    
	}


    return sc;
}


EventSelector::EventSelector( const std::string& name, ISvcLocator* svcloc ) 
: Service( name, svcloc)
, m_evtEnd(this,-1, 0) 
{
    declareProperty( "JobInput", m_jobInput=" ");
    m_inputDataList = new ListName; 
    // GlastEvtIterator m_evtEnd(this, -1, m_inputDataList->end());
}


EventSelector::~EventSelector() {
    delete m_inputDataList;
}

// IEvtSelector::setCriteria
// Create list of input files/Job Id's
StatusCode EventSelector::setCriteria( const std::string& criteria ) {
    
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
    rest      = rest.substr(lpos, -1);
    if( type == "IRFFILE" ) {
        m_criteriaType = IRFFILE;
    } else {
        log << MSG::ERROR << "Invalid Event Selection Criteria: " << criteria << endreq;
        return StatusCode::FAILURE;
    }
    sc = parseStringInList(rest, m_inputDataList);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Invalid Event Selection Criteria: " << criteria << endreq;
        return sc;
    }
    sc = serviceLocator()->getService ("GlastDetSvc",
         IID_IGlastDetSvc, reinterpret_cast<IInterface*&>( m_detSvc ));
    
    if(sc.isFailure()){
            log << MSG::ERROR << "Unable start Glast detector service" << endreq;
        return sc;
    }

    return sc;
}


// Parse criteria string: Fill in the list of input files or Job Id's
StatusCode EventSelector::parseStringInList( const std::string& namelist, ListName* inputDataList ) {
    std::string rest = namelist;
    while(true) {
        int ipos = rest.find_first_not_of(" ,");
        if (ipos == -1 ) break;
        rest = rest.substr(ipos, -1);            // remove blanks before
        int lpos  = rest.find_first_of(" ,");    // locate next blank
        if (lpos == -1 ) {
            rest = rest.substr(0,lpos );
            inputDataList->push_back( rest);       // insert last item in list and
            break;                                 // break
        }
        inputDataList->push_back( rest.substr(0,lpos ));   // insert in list
        rest = rest.substr(lpos, -1);                      // get the rest
    }
    return StatusCode::SUCCESS;
}

// IEvtSelector::setCriteria
StatusCode EventSelector::setCriteria( const SelectionCriteria& criteria ) {
    return StatusCode::SUCCESS;
}



// IEvtSelector::begin()
// Event Selector Iterator begin
// Called by the ApplicationMgr::intialize() method
IEvtSelector::Iterator* EventSelector::begin() const {
    MsgStream log(messageService(), name());
    StatusCode sc;
    // Ian Note Here is the problem I need to step through right here
    GlastEvtIterator* it = new GlastEvtIterator(this, -1, m_inputDataList->begin());
    
    if( (it->m_inputDataIt) ==  m_inputDataList->end()) { 
        return it;                               // If no data available  return begin = end
    }
    log << MSG::DEBUG << " First input data set is " << *(it->m_inputDataIt) << endreq;
    it->m_recId = 0;                            // If data file is found Initialize record iterator
    
    std::string fName;
    sc = getFileName(&(it->m_inputDataIt), fName);
        
    // assumes GlastDetectors have been init. Open for read
    //it->m_inputFile = new ResponseFile (fName.c_str(), true);
//	fName = "D:\\Gaudi2\\GlastSvc\\v1\\src\\test\\muontest.irf";//Ian Note
    sc = m_detSvc->openIRF(fName);
    if(sc.isFailure()){
        *(it) = m_evtEnd; 
        log << MSG::INFO << "failed to open file " << fName << endreq;
        return it;
    }

    
    (*it)++;  // increment to the first event
    
    return it;
}

//  Find out the name of the file from list of files or Jobs.
StatusCode EventSelector::getFileName(ListName::const_iterator* inputIt, std::string& fName) const {
    
    MsgStream log(messageService(), name());
    if( m_criteriaType == IRFFILE){                 // If CRITERIA = FILE Get File name
        fName = **inputIt;
    } else {
        log << MSG::ERROR << " Wrong Event selector Criteria " << endreq;
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}


// IEvtSelector::end()
// Event Selector Iterator end() function
IEvtSelector::Iterator* EventSelector::end() const {
    IEvtSelector::Iterator* it = (IEvtSelector::Iterator*)(&m_evtEnd);
    return it; 
}



// IEvtSelector::next
// Event Selector iterator next(it)
IEvtSelector::Iterator& EventSelector::next(IEvtSelector::Iterator& it) const {
    MsgStream log(messageService(), name());
    GlastEvtIterator* irfIt = dynamic_cast<GlastEvtIterator*>(&it);
    (irfIt->m_recId)++;
    log << MSG::DEBUG << "Reading Event " <<  irfIt->m_recId << endreq;
    
    // causes data from the current event in the IRF file to be loaded into the GlastDetector
    // objects
    StatusCode sc = m_detSvc->readIRF();// Ian Note: This is where the events are read in.

    if (sc.isFailure()){
        log << MSG::ERROR << " Failed to get Event " << irfIt->m_recId << endreq;
        log << MSG::ERROR << "The job will STOP normally after reading : " 
            << --(irfIt->m_recId) << " Events" << endreq;
        *(irfIt) = m_evtEnd;
    }

    return *irfIt;
}
// IEvtSelector::previous
IEvtSelector::Iterator& EventSelector::previous(IEvtSelector::Iterator& it) const {
    MsgStream log(msgSvc(), name());
    log << MSG::FATAL << " EventSelector Iterator, operator -- not supported " << endreq;
    
    return it;
}


//  IEvtSelector::reference

/*! Called from ApplicationMgr_i::nextEvent()       
Create root address and assign address to data service
IOpaqueAddress* addr = **m_evtIterator;
*/
IOpaqueAddress* EventSelector::reference(const IEvtSelector::Iterator& it) const 
{
    MsgStream log(msgSvc(), name());
    
    // convert to our iterator
    const GlastEvtIterator* irfIt = dynamic_cast<const GlastEvtIterator*>(&it);
    
    // get the file name and record id
    const std::string& fName = *(irfIt->m_inputDataIt);
    int recId = irfIt->m_recId;
    recId++;
    
    ListName::const_iterator* inputDataIt = (ListName::const_iterator*)(&irfIt->m_inputDataIt);
    log << MSG::DEBUG << "EventSelector::reference: File name is " << fName << endreq;
    
    // finally create an opaque address to pass back
    IOpaqueAddress* addr = DummyAddressFactory.instantiate(CLID_Event, fName, "PASS", recId);
    log << MSG::DEBUG << "EventSelector::reference: Created address:" << (void*)addr << endreq;
    return addr;
}

// IProperty::setProperty
StatusCode EventSelector::setProperty(const Property& property) {
    return StatusCode::SUCCESS;
}

// IProperty::getProperty
StatusCode EventSelector::getProperty(Property* property) const {
    return StatusCode::SUCCESS;
}

// IInterface::queryInterface
StatusCode EventSelector::queryInterface(const IID& riid, void** ppvInterface)  {
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


