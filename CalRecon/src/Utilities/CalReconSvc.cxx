// for class definition
#include "ICalReconSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

// for implementation
#include "Event/TopLevel/EventModel.h"

/**   
* @class CalReconSvc
*
* Data and features shared by all CalRecon actors.
*
* $Header$
*/

class CalReconSvc : public Service, public virtual ICalReconSvc
{
public:
    CalReconSvc( const std::string & name, ISvcLocator * svcLocator) : Service(name,svcLocator) {};
    virtual ~CalReconSvc() {};

    virtual StatusCode initialize() ;
    StatusCode queryInterface( const IID &, void** ppvUnknown ) ;
        
    //! update what depends on event state
    void reviewEvent() ;
    
    // generic services
    StatusCode              getStatus()         const {return m_status;}
    IDataProviderSvc*       getEventSvc()       const {return m_eventSvc;}
    IGlastDetSvc*           getDetSvc()         const {return m_detSvc;}

    // cal parameters
    int                     getCalNLayers()     const {return m_calNLayers;}
    double                  getCalCsIWidth()    const {return m_calCsIWidth;}
    double                  getCalCsIHeight()   const {return m_calCsIHeight;}

    // cal event data
    Event::CalXtalRecCol*   getXtalRecs()             {return m_calXtalRecCol;}
    Event::CalClusterCol*   getClusters()             {return m_calClusterCol;}

private :
    
    // some pointers to services  
    StatusCode            m_status ;
    IMessageSvc*          m_messageSvc ;
    IDataProviderSvc*     m_eventSvc ;
    IGlastDetSvc*         m_detSvc; 
    
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

StatusCode CalReconSvc::queryInterface( const IID & riid, void** ppvIF )
{
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

StatusCode CalReconSvc::initialize()
{
    Service::initialize() ;
     
    m_status = StatusCode::SUCCESS ;
  
    //svcLocator->service("MessageSvc",m_messageSvc,true) ;
    service("MessageSvc",m_messageSvc,true) ;
    MsgStream log(m_messageSvc,"CalReconSvc::CalReconSvc");

    if ((m_status = service("EventDataSvc",m_eventSvc,true)).isFailure())
    { 
        throw GaudiException("Service [EventDataSvc] not found", name(), m_status);
    }

    if ((m_status = service("GlastDetSvc",m_detSvc, true)).isFailure())
    { 
        throw GaudiException("Service [GlastDetSvc] not found", name(), m_status);
    }

    // cal const info
    if (!m_detSvc->getNumericConstByName(std::string("CALnLayer"), &m_calNLayers)) 
    { 
        m_status = StatusCode::FAILURE;
        throw GaudiException("Constant CALnLayer not defined", name(), m_status);
    }
    else
    {
        // David Chamont: number of layers is also hardcoded in CalCluster !
        //  it will not hurt if I check the two values are consistent.
        //Event::CalCluster fake(0,Point()) ;
        //if ( fake.getEneLayer().size() != ((unsigned int)m_calNLayers ))
        // { log<<MSG::FATAL<<"Inconsistent number of layers"<<endreq ; m_status = StatusCode::FAILURE ; }
    }

    if (!m_detSvc->getNumericConstByName(std::string("CsIWidth"),&m_calCsIWidth))
    { 
        m_status = StatusCode::FAILURE;
        throw GaudiException("Constant CsIWidth not defined", name(), m_status);
    }

    if (!m_detSvc->getNumericConstByName(std::string("CsIHeight"),&m_calCsIHeight))
    { 
        m_status = StatusCode::FAILURE;
        throw GaudiException("Constant CsIHeight not defined", name(), m_status);
    }

    return m_status ;
}
 
void CalReconSvc::reviewEvent()
{
    MsgStream log(m_messageSvc,"CalReconSvc::reviewEvent") ;

    // Ensure CalRecon/Event directory in TDS
    // David Chamont: I tried to move this section in the constructor above,
    //  which is called during the initization of CalClustersAlg, but
    //  for some not investigated reason it does not work.
    DataObject * pnode = 0 ;
    if ((getEventSvc()->retrieveObject(EventModel::CalRecon::Event,pnode)).isFailure()
      && (getEventSvc()->registerObject(EventModel::CalRecon::Event,new DataObject)).isFailure())
    { 
        log<<MSG::ERROR<<"Could not create CalRecon directory in TDS"<<endreq ; m_status = StatusCode::FAILURE ; 
    } 

    // CalXtalRecCol
    m_calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(getEventSvc(),EventModel::CalRecon::CalXtalRecCol) ; 
    if (!m_calXtalRecCol)
    { 
        log<<MSG::VERBOSE<<"No CalXtalRecCol"<<endreq ; 
    }
        
    // CalClusterCol
    m_calClusterCol = SmartDataPtr<Event::CalClusterCol>(getEventSvc(),EventModel::CalRecon::CalClusterCol) ;
    if (!m_calClusterCol)
    {
        m_calClusterCol = new Event::CalClusterCol() ;
        if ((getEventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,m_calClusterCol)).isFailure())
        {
            log<<MSG::ERROR<<"Cannot register CalClusterCol"<<endreq ;
        m_status = StatusCode::FAILURE ;
        }
    }

    return;
}


