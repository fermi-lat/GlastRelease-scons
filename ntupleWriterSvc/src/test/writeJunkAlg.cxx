#define ntupleWriter_writeJunkAlg_CPP 

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"


#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/INTuple.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

//------------------------------------------------------------------------------
/*! \class writeJunkAlg
\brief  alg to control writing of ntuples

  */

class writeJunkAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    writeJunkAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    std::string m_tupleName;
    float m_float;

    INTupleWriterSvc *ntupleWriteSvc;

};

//------------------------------------------------------------------------------
static const AlgFactory<writeJunkAlg>  Factory;
const IAlgFactory& writeJunkAlgFactory = Factory;
//------------------------------------------------------------------------------
/// 
writeJunkAlg::writeJunkAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){

    declareProperty("tupleName",  m_tupleName="");

}

//------------------------------------------------------------------------------
/*! 
*/
StatusCode writeJunkAlg::initialize() {
    
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get a pointer to our ntupleWriterSvc
    sc = serviceLocator()->getService("ntupleWriterSvc", 
        IID_INTupleWriterSvc, reinterpret_cast<IInterface*&>( ntupleWriteSvc));

    if( sc.isFailure() ) {
        log << MSG::ERROR << "writeJunkAlg failed to get the ntupleWriterSvc" << endreq;
        return sc;
    }

    return sc;
}



//------------------------------------------------------------------------------
StatusCode writeJunkAlg::execute() {
     
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    m_float = 700.05;
    sc = ntupleWriteSvc->addItem(m_tupleName.c_str(), "MyFirstItem", m_float);

    return sc;
}


//------------------------------------------------------------------------------
StatusCode writeJunkAlg::finalize() {
    StatusCode  sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize writeJunkAlg " << endreq;
 
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------



