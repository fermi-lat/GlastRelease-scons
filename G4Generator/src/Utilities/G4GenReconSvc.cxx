// for class definition
#include "GaudiKernel/Service.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

// for implementation
#include "G4GenErrorRecord.h"
#include "G4GenException.h"
#include "IG4GenErrorSvc.h"

/**   
* @class G4GenErrorSvc
*
* Data and features shared by all G4GenRecon actors.
*
* $Header$
*/

class G4GenErrorSvc : virtual public IG4GenErrorSvc, public Service, public virtual IIncidentListener 
{
public:
    G4GenErrorSvc( const std::string &, ISvcLocator * ) ;
    virtual ~G4GenErrorSvc() ;

    virtual StatusCode initialize() ;
    virtual StatusCode finalize() ;
    StatusCode queryInterface( const IID &, void** ) ;
        
    //! G4Genled when an event begin and end
    void handle( const Incident & ) ;
    
    //! register errors
    StatusCode handleError
          ( const std::string & catcherName,
            const std::string & comment ) ;

private :    
    // some pointers to services  
    MsgStream * m_log ;
    IDataProviderSvc * m_eventSvc ;
        
    // for errors
    bool m_saveBadEvents ;
    bool m_testExceptions ;
    bool m_printExceptions ;
    Event::EventHeader * m_header ;
    std::vector<G4GenErrorRecord *> m_errors ;
    double m_lastTime;
} ;

//==========================================================
// implementation
//==========================================================


DECLARE_SERVICE_FACTORY(G4GenErrorSvc) ;

StatusCode G4GenErrorSvc::queryInterface( const IID & riid, void** ppvIF ) 
{
    if (IID_IG4GenErrorSvc == riid)
    {
        *ppvIF = dynamic_cast<IG4GenErrorSvc*>(this) ;
        return StatusCode::SUCCESS ;
    }
    else
    { 
        return Service::queryInterface(riid,ppvIF) ; 
    }

}

G4GenErrorSvc::G4GenErrorSvc( const std::string & name, ISvcLocator * svcLocator )
                            : Service(name,svcLocator) 
{
    // If true, try very hard not to crash
    declareProperty("saveBadEvents",m_saveBadEvents=true) ;
    // turn on some bad things to see if exception handling is working
    declareProperty("testExceptions",m_testExceptions=false) ;
    // force exceptions to be printed
    declareProperty("printExceptions",m_printExceptions= false) ;
}

G4GenErrorSvc::~G4GenErrorSvc() {}


StatusCode G4GenErrorSvc::initialize() 
{    
    Service::initialize() ;
    setProperties() ; // said to be necessary for services
     
    StatusCode sc = StatusCode::SUCCESS ;
  
    //svcLocator->service("MessageSvc",m_messageSvc,true) ;
    IMessageSvc * messageSvc ;
    service("MessageSvc",messageSvc,true) ;
    m_log = new MsgStream(messageSvc,name()) ;

    if ((sc = service("EventDataSvc",m_eventSvc,true)).isFailure())
    { 
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    // incident listening
    IIncidentSvc * incsvc ;
    if ((sc = service("IncidentSvc",incsvc,true)).isFailure()) 
    { 
        throw GaudiException("Service [IncidentSvc] not found",name(),sc) ;
    }

    incsvc->addListener(this,"BeginEvent") ;
    incsvc->addListener(this,"EndEvent") ;

    // for errors
    m_lastTime = 0.0 ;

   return sc ;
   
}
 
void G4GenErrorSvc::handle( const Incident & inc ) 
{
    (*m_log) << MSG::DEBUG << "Incident " << inc.type() << endreq ;
    
    if (inc.type()=="BeginEvent") 
    {    
        // EventHeader
        m_header = SmartDataPtr<Event::EventHeader>(m_eventSvc,EventModel::EventHeader);
        if (!m_header) 
        {
            (*m_log) << MSG::ERROR << "Event header not found !" << endreq ;
        }

    } 
    else if (inc.type()=="EndEvent") 
    {
        m_lastTime = m_header->time();
    }
}

StatusCode G4GenErrorSvc::handleError(const std::string & catcherName, const std::string & comment )  
{        
    int run = m_header->run();
    int event  = m_header->event();
    G4GenErrorRecord * error = new G4GenErrorRecord( run,event,m_lastTime,catcherName,comment );
    m_errors.push_back(error);

    if (!m_printExceptions) error->print((*m_log) << MSG::DEBUG) ;
    else  error->print((*m_log) << MSG::WARNING) ;
    
    if (m_saveBadEvents) return StatusCode::SUCCESS;
    else return StatusCode::FAILURE ;
}

StatusCode G4GenErrorSvc::finalize() 
{    
    Service::finalize() ;
    
    int errorCount = m_errors.size() ;
    if (m_saveBadEvents&&((*m_log)<<MSG::INFO).isActive()) 
    {
        (*m_log) << endreq
                 << "====>> " << errorCount << (errorCount==1?" event":" events")
                 << " failed in this run" << endreq << endreq ;
    }

    std::vector<G4GenErrorRecord*>::iterator error ;
    for ( error=m_errors.begin() ; error!=m_errors.end() ; ++error ) 
    {
        (*error)->print((*m_log) << MSG::INFO);
    }

    delete m_log ; m_log = 0 ;
    return StatusCode::SUCCESS ;
}
