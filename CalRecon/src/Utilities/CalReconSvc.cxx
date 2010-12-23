
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
        StatusCode queryInterface( const InterfaceID &, void** ) ;
        
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
        double                  getCalCsILength()   const {return m_calCsILength;}
        double                  getCaltowerPitch()  const {return m_caltowerPitch;}
        bool                    getCalFlightGeom()  const {return m_calflightgeom;}

        // MomentsClusterInfo parameters.
        double getMciZeroSupprEnergy()              const { return m_mciZeroSupprEnergy; }
        double getMciXtalsTruncFrac()               const { return m_mciXtalsTruncFrac; }
        double getMciEneMomTruncFrac()              const { return m_mciEneMomTruncFrac; }

        // Quantities for the moments analysis.
        double getMaTransScaleFactor()              const { return m_maTransScaleFactor; }
        double getMaTransScaleFactorBoost()         const { return m_maTransScaleFactorBoost; }
        double getMaCoreRadius()                    const { return m_maCoreRadius; }
    
        // cal event data
        Event::CalXtalRecCol * getXtalRecs() ;
        Event::CalDigiCol *getDigis() ;
        Event::CalClusterCol * getClusters() ;

    private :
    
        // some pointers to services  
        MsgStream * m_log ;
        IDataProviderSvc * m_eventSvc ;
        IGlastDetSvc * m_detSvc; 
        
        // Quantities for the MomentsClusterInfo class.
        double m_mciZeroSupprEnergy;
        double m_mciXtalsTruncFrac;
        double m_mciEneMomTruncFrac;

        // Quantities for the moments analysis.
        double m_maTransScaleFactor;
        double m_maTransScaleFactorBoost;
        double m_maCoreRadius;

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
        //! CAL crystal length
        double m_calCsILength;
        //! CAL tower pitch
        double m_caltowerPitch;
        // !CAL tower numbers and flight mode
        int m_calnumX;
        int m_calnumY;
        bool m_calflightgeom;

        //! reconstructed data for crystals
        Event::CalXtalRecCol* m_calXtalRecCol ;
        //! the digis list
        Event::CalDigiCol* m_calDigiCol ;
        //! the clusters list
        Event::CalClusterCol* m_calClusterCol ;
} ;


//==========================================================
// implementation
//==========================================================


DECLARE_SERVICE_FACTORY(CalReconSvc) ;

StatusCode CalReconSvc::queryInterface( const InterfaceID & riid, void** ppvIF ) {

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

    // Energy correction (in MeV) for the zero suppression in the CAL.
    // This number is multiplied by the truncated number of xtals.
    declareProperty("MomentsClusterInfo_zeroSupprEnergy", m_mciZeroSupprEnergy=0.2) ;
    // Energy threshold (in fraction of the raw energy) used to count
    // the truncated number of xtals.
    declareProperty("MomentsClusterInfo_xtalsTruncFrac", m_mciXtalsTruncFrac=0.01) ;
    // Energy threshold (in fraction of the raw energy) used to calculate
    // the moments of the xtal energy distribution.
    declareProperty("MomentsClusterInfo_eneMomTruncFrac", m_mciEneMomTruncFrac=0.01) ;

    // Multiplicative factor driving the iterative moments analysis
    // (at each iteration the xtals whose distance exceeds this factor times the
    // transverse rms at the same step are dropped).
    declareProperty("CalMomentsAnalysis_transScaleFactor", m_maTransScaleFactor=1.0) ;
    // Multiplicative boost factor driving the change of the previous one
    // (transScaleFactor) at each iteration.
    declareProperty("CalMomentsAnalysis_transScaleFactorBoost", m_maTransScaleFactorBoost=2.0) ;
    // The radius (in terms of Moliere radii) of the cylinder around the axis
    // used to integrate the fraction of shower core energy in the moments analysis.
    declareProperty("CalMomentsAnalysis_coreRadius", m_maCoreRadius = 0.75) ;
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

    if (!m_detSvc->getNumericConstByName(std::string("CsILength"),&m_calCsILength))
    { 
        sc = StatusCode::FAILURE;
        throw GaudiException("Constant CsILength not defined", name(), sc);
    }

    if (!m_detSvc->getNumericConstByName(std::string("towerPitch"),&m_caltowerPitch))
    { 
        sc = StatusCode::FAILURE;
        throw GaudiException("Constant towerPitch not defined", name(), sc);
    }

    if (!m_detSvc->getNumericConstByName(std::string("xNum"),&m_calnumX))
    { 
        sc = StatusCode::FAILURE;
        throw GaudiException("Constant xNum not defined", name(), sc);
    }

    if (!m_detSvc->getNumericConstByName(std::string("yNum"),&m_calnumY))
    { 
        sc = StatusCode::FAILURE;
        throw GaudiException("Constant yNum not defined", name(), sc);
    }

    m_calflightgeom = true;
    if(m_calnumX==4 && m_calnumY==1) m_calflightgeom = false;
    
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
        m_calDigiCol = 0 ;
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

Event::CalDigiCol * CalReconSvc::getDigis() {
    
    if (m_calDigiCol!=0) return m_calDigiCol ;
    
    m_calDigiCol = SmartDataPtr<Event::CalDigiCol>(getEventSvc(),EventModel::Digi::CalDigiCol);
    if (!m_calDigiCol) { 
        (*m_log)<<MSG::VERBOSE<<"No CalDigiCol"<<endreq ; 
    }
        
    return m_calDigiCol ;
}

Event::CalClusterCol * CalReconSvc::getClusters() {
    
    if (m_calClusterCol!=0) return m_calClusterCol ;
    
    // Ensure CalRecon/Event directory in TDS
    //DataObject * pnode = 0 ;
    //if ((getEventSvc()->retrieveObject(EventModel::CalRecon::Event,pnode)).isFailure()
    //  && (getEventSvc()->registerObject(EventModel::CalRecon::Event,new DataObject)).isFailure()) { 
    //    throw CalException("cannot register Event/CalRecon") ;
    //} 

    // CalClusterCol
    m_calClusterCol = SmartDataPtr<Event::CalClusterCol>(getEventSvc(),EventModel::CalRecon::CalClusterCol) ;
    //if (!m_calClusterCol) {
    //    m_calClusterCol = new Event::CalClusterCol() ;
    //    if ((getEventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,m_calClusterCol)).isFailure()) {
    //        throw CalException("cannot register CalClusterCol") ;
    //    }
    //}
    
    return m_calClusterCol ;

}


