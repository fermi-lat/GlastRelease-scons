/** @file Count.cxx

*/
#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataSvc.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"


/** @class Count
 * @brief Alg that counts and records the value in the jobinfo tuple
 * 
 * $Header$
*/
class Count : public Algorithm
{   
public:
    Count(const std::string& name, ISvcLocator* pSvcLocator); 
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
private:
    INTupleWriterSvc* m_ntupleWriterSvc;
    int m_count;
    
    
};

//  factory stuff
static const AlgFactory<Count>  Factory;
const IAlgFactory& CountFactory = Factory;


Count::Count(const std::string& name, ISvcLocator* pSvcLocator) 
:Algorithm(name, pSvcLocator) 
, m_count(0)
{
}


StatusCode Count::initialize() { 
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    StatusCode sc =StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    sc = service("RootTupleSvc", m_ntupleWriterSvc);
    
    if(sc.isFailure())    {
        log << MSG::ERROR << "Could not locate the ntupleSvc" <<endreq;
    }

    return sc;
    
}


StatusCode Count::execute() {

    ++m_count;
    MsgStream   log( msgSvc(), name() );
    StatusCode  sc = StatusCode::SUCCESS;

    return sc;
}


StatusCode Count::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "counted " << m_count << " times" << endreq;
    m_ntupleWriterSvc->addItem("jobinfo", name(), &m_count);
    return StatusCode::SUCCESS;
}

