#define AcdDigi_ReadDigi_CPP 

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "GaudiKernel/ObjectVector.h"

#include "GaudiKernel/Algorithm.h"

/** @class ReadDigi.cpp
 * @brief AcdDigi Test algorithm
 * $Header$
 */

class ReadDigi : public Algorithm {

public:
  ReadDigi(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:

};

static const AlgFactory<ReadDigi>  Factory;
const IAlgFactory& ReadDigiFactory = Factory;

ReadDigi::ReadDigi(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {

}

StatusCode ReadDigi::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    return StatusCode::SUCCESS;
}


StatusCode ReadDigi::execute() {
    using namespace Event;
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "execute" << endreq;


    return sc;
}


StatusCode ReadDigi::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    
    return StatusCode::SUCCESS;
}






