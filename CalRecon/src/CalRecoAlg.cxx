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

#include "GlastEvent/data/CalData.h"

#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
#define TUPLE 1

static const AlgFactory<CalRecoAlg>  Factory;
const IAlgFactory& CalRecoAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.
CalRecoAlg::CalRecoAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator), m_detSvc(0)
, m_recon(0), m_gsummary(0) {
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
    StatusCode sc = service("GlastDetSvc", m_detSvc);
    
    
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
        return sc;
    }
    
    // test: get a constant from the ini file
    m_ini = const_cast<xml::IFile*>(m_detSvc->iniFile()); //OOPS!
    int nx = m_ini->getInt("glast", "xNum");
    
    m_recon=new CalRecon;
    
    // get the Gui service (not required)
    IGuiSvc* guiSvc=0;
    sc = service("GuiSvc", guiSvc);

    if (!sc.isSuccess ()){
        log << MSG::WARNING << "No GuiSvc, so no display" << endreq;
        sc =StatusCode::SUCCESS; 
    }else
    {
       guiSvc->guiMgr()->display().add(m_recon->displayRep(), "Cal reco");
    }

#ifdef TUPLE    
    // create and set the entries into a Gaudi Tuple
    m_gsummary = new SummaryData<GaudiGlastTuple>(
        *new GaudiGlastTuple(ntupleSvc(),"Cal Recon output", "CALRECON"));
    m_recon->accept(*m_gsummary);           
#endif
    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalRecoAlg::execute() {
    
 
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "execute" << endreq;
    
    // get the CsiData object from the TDS by a converter
    //static std::string path("/Event/Raw/TdCsIDatas");
    static std::string path("/Event/data/CalData");
    SmartDataPtr<CalData> csi(eventSvc(), path);
    
    if( 0==csi) { log << MSG::ERROR << "could not find \""<< path <<"\"" << endreq;
       return StatusCode::FAILURE;
    }
        
    // see what is there
    //csi->printOn(std::cout);
    
    // create the CalRecon object from the reconstrution package and pass data to it.
    
    m_recon->clear();
    m_recon->reconstruct(csi);
    
    // print out the  tuple
    //m_recon->accept(PrintReconData(std::cout));
    
#ifdef TUPLE
    // fill and write out the Gaudi tuple   
    m_gsummary->tuple()->fill();

#endif
    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalRecoAlg::finalize() {
    StatusCode  sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    delete m_recon;
    delete m_gsummary;

    //try to read tuple back
    NTuplePtr nt(ntupleSvc(), "/NTUPLES/CALRECON/1");
    if( nt) {
        int count=0;
        NTuple::Item<float>esum;
        sc = nt->item("CsI_Energy_Sum", esum);
        while( ntupleSvc()->readRecord(nt.ptr()).isSuccess() ) {
            log << MSG::INFO << " Entry [" << (count++) <<"]" << "Etot=" << esum << endreq;
        }
        log << MSG::INFO << count << " entries " << endreq;
               
    }else {
        log << MSG::INFO << "Could not open tuple  /NTUPLES/CalRecon/1" << endreq;
    }
    
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
/** Can be used to check to see if the newstyle NTuple was written properly.
*/
StatusCode CalRecoAlg::printNewNTuple(std::string dname) {
    StatusCode status;
    MsgStream log(msgSvc(), name());
    
    NTuplePtr nt = m_gsummary->tuple()->getNTuple();
    if( nt) {
        NTuple::Item<float> data;
        status = nt->item(dname ,data);
       
        if(status.isSuccess()) {
            log << MSG::INFO << "Test Value of :"
                << dname << "= " << data << endreq;
        } else {
            log << MSG::ERROR << "Didn't load the Item!" << endreq;
        }
    } else {
        log << MSG::ERROR << "Dead NTuple !" << endreq;
    }
    return status;
}





