/** @file OrbitSvc.cxx
@brief declaration and definition of the class OrbitSvc

$Header$

*/

// Include files
// Gaudi system includes
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/Property.h"


#include "facilities/Util.h"
#include "facilities/Timestamp.h"

#include "astro/GPS.h"


/** 
* \class OrbitSvc
*
*/

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IOrbitSvc("OrbitSvc", 0, 0); 


class OrbitSvc : virtual public Service  {
public:
    OrbitSvc(const std::string& name, ISvcLocator* pSvcLocator);


    //------------------------------------------------------------------
    //  stuff required by a Service

    /// perform initializations for this service. 
    virtual StatusCode initialize ();

    /// perform the finalization, as required for a service.
    virtual StatusCode finalize ();

    /// Query interface
    virtual StatusCode queryInterface( const InterfaceID& riid, void** ppvUnknown );

private: 

    StringArrayProperty m_pointingHistory;///< history file name and launch date
    DoubleArrayProperty m_pointingDirection; ///< (ra, dec) for pointing
    DoubleProperty m_zenithTheta; ///< set for zenith
    DoubleProperty m_rocking_angle; // x-axis

};
//------------------------------------------------------------------------


//static const SvcFactory<OrbitSvc>  Factory;
//const ISvcFactory& OrbitSvcFactory = Factory;
DECLARE_SERVICE_FACTORY(OrbitSvc);

//------------------------------------------------------------------------
//! ctor
OrbitSvc::OrbitSvc(const std::string& name, ISvcLocator* svc)
: Service(name,svc)
{

    declareProperty("PointingHistory",   m_pointingHistory); // doublet, filename and launch date
    declareProperty("pointingDirection", m_pointingDirection);
    declareProperty("zenithTheta",       m_zenithTheta=-99);
    declareProperty("rocking_angle", m_rocking_angle=0); // set non-zero to enable rocking



}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode OrbitSvc::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    using astro::GPS;

    // Use the Job options service to set the Algorithm's parameters
    setProperties();


    //set the input file to be used as the pointing database, if used
    if(! m_pointingHistory.value().empty()){
        std::string filename(m_pointingHistory.value()[0]);
        facilities::Util::expandEnvVar(&filename);
        double offset = 0;
        bool horizontalflag(false);
        if( m_pointingHistory.value().size()>1){
            std::string field(m_pointingHistory.value()[1]);
            if(! field.empty() ) { // allow null string
                facilities::Timestamp jt(m_pointingHistory.value()[1]);
                offset = (astro::JulianDate(jt.getJulian())-astro::JulianDate::missionStart())*astro::JulianDate::secondsPerDay;
            }
        }

        if( m_pointingHistory.value().size()>2){
            std::string field(m_pointingHistory.value()[2]);
            horizontalflag =! field.empty();
        }
        log << MSG::INFO << "Loading Pointing History File : " << filename 
            << " with MET offset "<< offset <<  endreq;
        if( horizontalflag){
            log << MSG::INFO << "Will override x-direction to be horizontal"<<endreq;
        }

        astro::GPS::instance()->setPointingHistoryFile(filename, offset, horizontalflag);
    }
    else if( m_pointingDirection.value().size()==2 ) {
        // if direction set, ignore everything else!
        double ra(m_pointingDirection.value()[0]), dec(m_pointingDirection.value()[1]);
        astro::GPS::instance()->setPointingDirection(astro::SkyDir(ra, dec));
        log << MSG::INFO << "set to point at ra,dec= " << ra << ", "<<dec << endreq;

    }else if( m_zenithTheta>=0) {
        // if set from defalt, set to this value
        //m_fluxSvc->setRockType(astro::GPS::EXPLICIT, m_zenithTheta);
        GPS::instance()->setRockType(GPS::EXPLICIT);
        GPS::instance()->rockingDegrees(m_zenithTheta);

        log << MSG::INFO << "set to zenith angle " << m_zenithTheta << " degrees" << endreq;
    }else if(m_rocking_angle >0 ){
        //output to record the pointing settings
        //then this line sets the rocking type, as well as the rocking angle.
        //m_fluxSvc->setRockType(GPS::ONEPERORBIT ,m_rocking_angle);
        GPS::instance()->setRockType(GPS::ONEPERORBIT);
        GPS::instance()->rockingDegrees(m_zenithTheta);
        log << MSG::INFO << "Once per orbit rocking Angle: " << m_rocking_angle << " degrees" << endreq;

    }


    return sc;
}

/// Query interface
StatusCode OrbitSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
    if ( IID_IOrbitSvc.versionMatch(riid) )  {
        *ppvInterface = (OrbitSvc*)this;
    } else  {
        return Service::queryInterface(riid, ppvInterface);
    }

    addRef();
    return SUCCESS;
}

StatusCode OrbitSvc::finalize()
{
    return StatusCode::SUCCESS;
}

