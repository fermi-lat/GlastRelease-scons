// $Header$

// Include files
#include "CalRecon/CalRecoAlg.h"

#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "xml/IFile.h"

#include "reconstruction/GlastTuple.h"
#include "reconstruction/PrintReconData.h"
#include "reconstruction/SummaryData.h"
#include "reconstruction/GlastTuple.h"

#include "GlastEvent/Raw/TdCsIData.h"

static const AlgFactory<CalRecoAlg>  Factory;
const IAlgFactory& CalRecoAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.
CalRecoAlg::CalRecoAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator), m_detSvc(0) {
}


//------------------------------------------------------------------------------
/*! The "functional" part of the class: For the EmptyAlgorithm example they do
nothing apart from print out info messages.
NB in the initialize method: you must explicitly initialize the base class
before using any services (message service, event data service etc.) otherwise 
the behaviour will be unpredictable at best.
*/
StatusCode CalRecoAlg::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // now try to find the GlastDevSvc service
#if 1  // old way, Gaudi v5.  
    IGlastDetSvc* detSvc = 0;
    const IID&  IID_IGlastDetSvc  =  401; // wired it for now!
    StatusCode sc = serviceLocator()->getService ("GlastDetSvc",
        IID_IGlastDetSvc, reinterpret_cast<IInterface*&>( detSvc ));
    
    m_detSvc = detSvc;

#else //this is the new way in v6, 
    StatusCode sc = service("GlastDetSvc", m_detSvc);
#endif    

    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
    }

    // test: get a constant from the ini file
    m_ini = const_cast<xml::IFile*>(m_detSvc->iniFile()); //OOPS!
    int nx = m_ini->getInt("glast", "xNum");

    m_recon=new CalRecon;

    // define the tuple
    m_summary = new  SummaryData<GlastTuple>(*new GlastTuple("test cal tuple")) ;
    m_recon->accept(*m_summary);

    m_summary->tuple()->writeHeader(std::cout);
    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalRecoAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "execute" << endreq;
        
    // get the CsiData object from the TDS by a converter
    SmartDataPtr<TdCsIData> csi(eventSvc(), "/Event/Raw/TdCsIDatas");

    // see what is there
    csi->printOn(std::cout);

    // create the CalRecon object from the reconstrution package and pass data to it.

    m_recon->clear();
    m_recon->reconstruct(csi);

    // print out the  tuple
    m_recon->accept(PrintReconData(std::cout));

    // fill the tuple and print the line
    m_summary->tuple()->fill();
    std::cout << *(m_summary->tuple());

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalRecoAlg::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    delete m_recon;
    delete m_summary;
    
    return StatusCode::SUCCESS;
}






