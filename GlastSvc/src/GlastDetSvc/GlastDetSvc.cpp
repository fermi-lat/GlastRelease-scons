// $Header$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#include "GlastDetSvc.h"

#include <fstream>
#include "Gaudi/Kernel/SvcFactory.h"
#include "Gaudi/MessageSvc/MsgStream.h"
#include "instrument/Instrument.h"
#include "xml/IFile.h"

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
, m_irfFile ("")
{
    
    // declare the properties
    declareProperty ("PersistencyFile", m_xmlFile);
    declareProperty ("IniFile", m_iniFile);
    declareProperty ("IRFFile", m_irfFile);
    
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
#if 0
    else if (IID_IGlastIRFLoadSvc == riid) {
        *ppvIF = dynamic_cast<IGlastIRFLoadSvc*> (this);
        return StatusCode::SUCCESS;
    }
#endif
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
    
    // now create and initialize the instruemtn
    m_instrument = new Instrument;
    if (m_instrument->initialize(m_iniFile, m_xmlFile) ) status=StatusCode::FAILURE;
    
    return status;
}

// finalize
StatusCode GlastDetSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    
    // get rid of a glast detector hierarchy
    delete m_instrument;
    m_instrument = 0;
    
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
    
    int rc = m_instrument->openIRF(m_irfFile);
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
extern void __cdecl FATAL ( const char* msg )
{
    (*s_log) << MSG::FATAL <<  msg << endreq;
}

//! Satisfy external
extern void __cdecl WARNING ( const char* msg )
{
    (*s_log) << MSG::WARNING << msg << endreq;
}

