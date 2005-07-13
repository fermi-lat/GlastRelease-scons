
// for class definition
#include <CalRecon/ICalReconSvc.h>
#include "Event/TopLevel/Event.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"

// for implementation
#include "CalErrorRecord.h"
#include "Event/TopLevel/EventModel.h"
#include "src/Utilities/CalException.h"

/**   
* @class CalReconSvc
*
* Data and features shared by all CalRecon actors.
*
* $Header$
*/

class CalReconSvc
  : public Service,
    public virtual ICalReconSvc,
    public virtual IIncidentListener {

    public:

        CalReconSvc( const std::string &, ISvcLocator * ) ;
        virtual ~CalReconSvc() ;

        virtual StatusCode initialize() ;
        virtual StatusCode finalize() ;
        StatusCode queryInterface( const IID &, void** ) ;
        
        //! called when an event begin and end
        void handle( const Incident & ) ;
    
        //! register errors
        StatusCode handleError
          ( const std::string & catcherName,
            const std::string & comment ) ;

        // generic services
        IDataProviderSvc*       getEventSvc()       const {return m_eventSvc;}
        IGlastDetSvc*           getDetSvc()         const {return m_detSvc;}

        // cal parameters
        int                     getCalNLayers()     const {return m_calNLayers;}
        double                  getCalCsIWidth()    const {return m_calCsIWidth;}
        double                  getCalCsIHeight()   const {return m_calCsIHeight;}
    
        // cal event data
        Event::CalXtalRecCol * getXtalRecs() ;
        Event::CalClusterCol * getClusters() ;

    private :
    
        // some pointers to services  
        MsgStream * m_log ;
        IDataProviderSvc * m_eventSvc ;
        IGlastDetSvc * m_detSvc; 
        
        // for errors
        bool m_saveBadEvents ;
        bool m_testExceptions ;
        bool m_printExceptions ;
        Event::EventHeader * m_header ;
        std::vector<CalErrorRecord *> m_errors ;
        double m_lastTime;
    
        //! CAL number of layers
        int                   m_calNLayers ;
        //! CAL crystal width
        double                m_calCsIWidth ;
        //! CAL crystal height
        double                m_calCsIHeight ;
        
        //! reconstructed data for crystals
        Event::CalXtalRecCol* m_calXtalRecCol ;
        //! the clusters list
        Event::CalClusterCol* m_calClusterCol ;
} ;


//==========================================================
// implementation
//==========================================================


DECLARE_SERVICE_FACTORY(CalReconSvc) ;

StatusCode CalReconSvc::queryInterface( const IID & riid, void** ppvIF ) {

    if (IID_ICalReconSvc == riid)
    {
        *ppvIF = dynamic_cast<ICalReconSvc*>(this) ;
        return StatusCode::SUCCESS ;
    }
    else
    { 
        return Service::queryInterface(riid,ppvIF) ; 
    }

}

CalReconSvc::CalReconSvc( const std::string & name, ISvcLocator * svcLocator )
  : Service(name,svcLocator) {

    // If true, try very hard not to crash
    declareProperty("saveBadEvents",m_saveBadEvents=true) ;
    // turn on some bad things to see if exception handling is working
    declareProperty("testExceptions",m_testExceptions=false) ;
    // force exceptions to be printed
    declareProperty("printExceptions",m_printExceptions= false) ;

}

CalReconSvc::~CalReconSvc() {
}


StatusCode CalReconSvc::initialize() {
    
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

    if ((sc = service("GlastDetSvc",m_detSvc, true)).isFailure())
    { 
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }

    // incident listening
    IIncidentSvc * incsvc ;
    if ((sc = service("IncidentSvc",incsvc,true)).isFailure()) { 
        throw GaudiException("Service [IncidentSvc] not found",name(),sc) ;
    }
    incsvc->addListener(this,"BeginEvent") ;
    incsvc->addListener(this,"EndEvent") ;

    // cal const info
    if (!m_detSvc->getNumericConstByName(std::string("CALnLayer"), &m_calNLayers)) 
    { 
        sc = StatusCode::FAILURE;
        throw GaudiException("Constant CALnLayer not defined", name(), sc);
    }
    else
    {
        // David Chamont: number of layers is also hardcoded in CalCluster !
        //  it will not hurt if I check the two values are consistent.
        //Event::CalCluster fake(0,Point()) ;
        //if ( fake.getEneLayer().size() != ((unsigned int)m_calNLayers ))
        // { (*m_log)<<MSG::FATAL<<"Inconsistent number of layers"<<endreq ; sc = StatusCode::FAILURE ; }
    }

    if (!m_detSvc->getNumericConstByName(std::string("CsIWidth"),&m_calCsIWidth))
    { 
        sc = StatusCode::FAILURE;
        throw GaudiException("Constant CsIWidth not defined", name(), sc);
    }

    if (!m_detSvc->getNumericConstByName(std::string("CsIHeight"),&m_calCsIHeight))
    { 
        sc = StatusCode::FAILURE;
        throw GaudiException("Constant CsIHeight not defined", name(), sc);
    }

   // for errors
   m_lastTime   = 0.0 ;

   return sc ;
   
}
 
void CalReconSvc::handle( const Incident & inc ) {

    (*m_log)<<MSG::DEBUG<<"Incident "<<inc.type()<<endreq ;
    
    if (inc.type()=="BeginEvent") {
        
        // EventHeader
        m_header = SmartDataPtr<Event::EventHeader>(getEventSvc(),EventModel::EventHeader) ;
        if (!m_header) {
            (*m_log)<<MSG::ERROR<<"Event header not found !"<<endreq ;
        }
        
        // reset other values, which cannot be prepared as long as
        // reading input files is done with an algo.
        m_calXtalRecCol = 0 ;
        m_calClusterCol = 0 ;

    } else if (inc.type()=="EndEvent") {

        m_lastTime = m_header->time() ;

    }

}

StatusCode CalReconSvc::handleError
  ( const std::string & catcherName,
    const std::string & comment )  {
        
    int run = m_header->run() ;
    int event  = m_header->event() ;
    CalErrorRecord * error = new CalErrorRecord
      ( run,event,m_lastTime,catcherName,comment ) ;
    m_errors.push_back(error) ;

    if (!m_printExceptions) error->print((*m_log)<<MSG::DEBUG) ;
    else  error->print((*m_log)<<MSG::WARNING) ;
    
    if (m_saveBadEvents) return StatusCode::SUCCESS;
    else return StatusCode::FAILURE ;

}

StatusCode CalReconSvc::finalize() {
    
    Service::finalize() ;
    
    int errorCount = m_errors.size() ;
    if (m_saveBadEvents&&((*m_log)<<MSG::INFO).isActive()) {
        (*m_log)<<endreq
          <<"====>> "<<errorCount<<(errorCount==1?" event":" events")
          <<" failed in this run"<<endreq<<endreq ;
    }
    std::vector<CalErrorRecord*>::iterator error ;
    for ( error=m_errors.begin() ; error!=m_errors.end() ; ++error ) {
        (*error)->print((*m_log)<<MSG::INFO) ;
    }

    delete m_log ; m_log = 0 ;
    return StatusCode::SUCCESS ;
}
 
Event::CalXtalRecCol * CalReconSvc::getXtalRecs() {
    
    if (m_calXtalRecCol!=0) return m_calXtalRecCol ;
    
    m_calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(getEventSvc(),EventModel::CalRecon::CalXtalRecCol) ; 
    if (!m_calXtalRecCol) { 
        (*m_log)<<MSG::VERBOSE<<"No CalXtalRecCol"<<endreq ; 
    }
        
    return m_calXtalRecCol ;
}

Event::CalClusterCol * CalReconSvc::getClusters() {
    
    if (m_calClusterCol!=0) return m_calClusterCol ;
    
    // Ensure CalRecon/Event directory in TDS
    DataObject * pnode = 0 ;
    if ((getEventSvc()->retrieveObject(EventModel::CalRecon::Event,pnode)).isFailure()
      && (getEventSvc()->registerObject(EventModel::CalRecon::Event,new DataObject)).isFailure()) { 
        throw CalException("cannot register Event/CalRecon") ;
    } 

    // CalClusterCol
    m_calClusterCol = SmartDataPtr<Event::CalClusterCol>(getEventSvc(),EventModel::CalRecon::CalClusterCol) ;
    if (!m_calClusterCol) {
        m_calClusterCol = new Event::CalClusterCol() ;
        if ((getEventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,m_calClusterCol)).isFailure()) {
            throw CalException("cannot register CalClusterCol") ;
        }
    }
    
    return m_calClusterCol ;

}


