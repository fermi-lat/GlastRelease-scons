// Gaudi 
#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataSvc.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"



//------------------------------------------------------------------------------
/*! \class WriteTupleAlg
Save the status of the HepRandomEngine used.
<br>

*/
class WriteTupleAlg : public Algorithm
{   
public:
    WriteTupleAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
private:
    INTupleWriterSvc* m_ntupleWriterSvc;
    int m_count;
    
    
};

//------------------------------------------------------------------------------
//  factory stuff
static const AlgFactory<WriteTupleAlg>  Factory;
const IAlgFactory& WriteTupleAlgFactory = Factory;

//------------------------------------------------------------------------------

WriteTupleAlg::WriteTupleAlg(const std::string& name, ISvcLocator* pSvcLocator) 
:Algorithm(name, pSvcLocator) 
{
    m_count =0;
}


StatusCode WriteTupleAlg::initialize() { 
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    StatusCode sc =StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    

    sc = service("ntupleWriterSvc", m_ntupleWriterSvc);
    
    if(sc.isFailure())
    {
        log << MSG::ERROR << "Could not locate the ntupleSvc" <<endreq;
    }


    return sc;
    
}


StatusCode WriteTupleAlg::execute() {

    
    MsgStream   log( msgSvc(), name() );
    StatusCode  sc = StatusCode::SUCCESS;

    sc = m_ntupleWriterSvc->saveNTuples();
    if(sc.isFailure())
    {
        log << MSG::ERROR << "Failed to save the ntuples" << endreq;
    }
    /*
        m_count++;
    if(m_count == 10)
    {
        return StatusCode::FAILURE;
    }
    */

    return sc;
}


//------------------------------------------------------------------------------
StatusCode WriteTupleAlg::finalize() {
    
    MsgStream log(msgSvc(), name());
    return StatusCode::SUCCESS;
}

