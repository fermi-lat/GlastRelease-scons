// $Header$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#include "GlastDetSvc.h"

#include <fstream>
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "instrument/Instrument.h"
#include "xml/IFile.h"

#include "GlastSvc/GlastDetSvc/IGeometry.h"
#include "GlastSvc/GlastDetSvc/IMedia.h"

#include "DMmanager.h"
#include "MediaVisitor.h"
#include "GeometryVisitor.h"

// define the class-id (unique) for the GlastDetSvc
// moved to Kernel! const IID&  IID_IGlastDetSvc  =  401;   // Unique to GLAST 

// define the class-id (unique) for the GlastIRFSvc
//const IID&  IID_IGlastIRFLoadSvc  =  402;   // Unique to GLAST


// declare the service factories for the GlastDetSvc
static SvcFactory<GlastDetSvc> a_factory;
const ISvcFactory& GlastDetSvcFactory = a_factory;


// file scope for WARNING and ERROR
static MsgStream * s_log;

// external declarations of the service identifiers
// HSG - not necessary from the Gaudi example - extern const IID& IID_IGlastDetSvc;

// ------------------------------------------------
// Implementation of the GlastDetSvc class
// ------------------------------------------------
/// Standard Constructor
GlastDetSvc::GlastDetSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc)
, m_instrument (0)
, m_xmlFile ("")
, m_iniFile ("")
, m_dm(0)
{
    
    // declare the properties
    declareProperty ("PersistencyFile", m_xmlFile);
    declareProperty ("IniFile", m_iniFile);

    // for DetModel visitor
    declareProperty("xmlfile",     m_xmlfile="-");
    declareProperty("topVolume",   m_topvol="LAT");
    declareProperty("visitorMode", m_visitorMode="propagate");


    
}
/// Standard Destructor
GlastDetSvc::~GlastDetSvc()  
{
    // this is  repetetive in case finalize was
    // not called.
    
    finalize();
}

// queryInterface
StatusCode  GlastDetSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_IGlastDetSvc == riid) {
        *ppvIF = dynamic_cast<IGlastDetSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}


// initialize
StatusCode GlastDetSvc::initialize () 
{
    StatusCode  status = StatusCode::SUCCESS;
    
    // call super-class
    Service::initialize ();
    
    // bind all of the properties for this service
    status = setProperties ();
    
    
    // open the message log
    MsgStream log( msgSvc(), name() );
    s_log = & log;  // make available globally while executing the following
    
    log << MSG::DEBUG << "Loading instrument ...";
    // now create and initialize the instrurment
    m_instrument = new Instrument;
    if (m_instrument->initialize(m_iniFile, m_xmlFile) ) status=StatusCode::FAILURE;
    log << MSG::DEBUG << "done. Loaded "<< m_instrument->detector_count() << " detectors." << endreq;

    
    // setup the detModel geometry, so can be visited below
    m_dm = new DMmanager;
    m_dm->init("-",m_visitorMode,m_topvol);
    log << MSG::INFO;
    m_dm->printSetup(log.stream());
    log << endreq;

    return status;
}
void GlastDetSvc::setDetector(GlastDetector* d)
{ 
    m_instrument->setDetector(d);    
}

// finalize
StatusCode GlastDetSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    
    // get rid of a glast detector hierarchy
    delete m_instrument;
    m_instrument = 0;
    
    // remove the detModel interface.
    delete m_dm;
    m_dm=0;
    return status;
}


// loadNextEvent - attempt to load the next event from the
//                 IRF file.
StatusCode  GlastDetSvc::readIRF ()
{
    StatusCode  status = StatusCode::SUCCESS;
    
    if (m_instrument->readIRF() ) status = StatusCode::FAILURE;
    return status;
}


// access the type of this service
const IID&  GlastDetSvc::type () const {
    return IID_IGlastDetSvc;
}



// loadIRFFile - attempt to open the IRF file here
StatusCode  GlastDetSvc::openIRF (std::string filename)
{
    MsgStream log( msgSvc(), name() );
    s_log = & log;  // make available globally while executing the following
    
    int rc = m_instrument->openIRF(filename);
    return rc==0? StatusCode::SUCCESS : StatusCode::FAILURE;
}

//! return the number of detectors currently with data
int GlastDetSvc::detectors_with_data()const 
{
    return m_instrument->detectors_with_data();
}

GlastDetector* GlastDetSvc::getRootDetector(void)
{
    return m_instrument->rootDetector();
}

//! accept a visitor to traverse the structure note that must be on a const object
void GlastDetSvc::accept(DetectorConverter& v)const {
    m_instrument->accept(v);
}

//! access to the IFile containing detector constants (implements IGlastDetSvc)
const xml::IFile* GlastDetSvc::iniFile()const{
    return m_instrument->iniFile();
}

///! satisfy external. Need to make it really fatal!
extern void FATAL ( const char* msg )
{
    (*s_log) << MSG::FATAL <<  msg << endreq;
}

//! Satisfy external
extern void WARNING ( const char* msg )
{
    (*s_log) << MSG::WARNING << msg << endreq;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//            the new detModel interface calls
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetSvc::accept(IGeometry& geom)
{
    m_dm->accept(&GeometryVisitor(geom));
}

void GlastDetSvc::accept(IMedia& media)
{
    m_dm->accept(&MediaVisitor(media));
}
