// $Header$

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include <string>


/*! \class G4TestAlg
\brief 
  This is a place to put test code to examine what G4Generator has done.
  */

class G4TestAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    G4TestAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    
};


static const AlgFactory<G4TestAlg>  Factory;
const IAlgFactory& G4TestAlgFactory = Factory;

//
G4TestAlg::G4TestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){
 
}


/*! */
StatusCode G4TestAlg::initialize() {
      

    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initializing..." << endreq;
    StatusCode  sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    return sc;
}


StatusCode G4TestAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    return sc;
}


StatusCode G4TestAlg::finalize() {
    
    return StatusCode::SUCCESS;
}






