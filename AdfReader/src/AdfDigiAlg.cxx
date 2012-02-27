// Gaudi specific include files
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"


#include "AdfEvent/AdfEvent.h"
#include "AncillaryDataEvent/Digi.h"

const std::string TDSdir = "/Event/AncillaryEvent";

class AdfDigiAlg : public Algorithm
{

public:
  AdfDigiAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~AdfDigiAlg() {}

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  StatusCode RegisterDigi();

private:
    IDataProviderSvc    *m_dataSvc;
  
};

//static const AlgFactory<AdfDigiAlg>  Factory;
//const IAlgFactory& AdfDigiAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(AdfDigiAlg);

AdfDigiAlg::AdfDigiAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : Algorithm(name, pSvcLocator) 
{
}

StatusCode AdfDigiAlg::initialize()
{
 // Purpose and Method:  Called once before the run begins.
  // This method:
  // 1) locates the Event Data Server.
  // 2) creates the AncillaryDataServer
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << "initialize" << endreq;
  setProperties();
  IService* iService = 0;
  sc = serviceLocator()->service("EventDataSvc", iService, true);
  
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "could not find DataSvc !" << endreq;
    return sc;
  }
  m_dataSvc = dynamic_cast<IDataProviderSvc*>(iService);
  return sc;  
}


StatusCode AdfDigiAlg::execute()
{  
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << "execute" << endreq;

  sc = RegisterDigi();
  if (sc.isFailure()) {
    log << MSG::INFO << "Failed to register or find AdfDigi" << endreq;
    return StatusCode::SUCCESS;
  }
  return sc;
}


StatusCode AdfDigiAlg::RegisterDigi()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  // Get the ADF event:
  std::string TDSobj = TDSdir + "/AdfEvent";
  log << MSG::DEBUG << " Get "<<TDSobj<<" from the TDS " << endreq;

  SmartDataPtr<AncillaryData::AdfEvent> currentEvent(eventSvc(),TDSobj);
  if( !currentEvent) 
    {
      log<<MSG::INFO<<"No AdfEvent found on TDS !" <<endreq;
      return sc;
    }
  
  // CREATE A DIGI EVENT
  AncillaryData::Digi *digiEvent = new AncillaryData::Digi(currentEvent);
  //digiEvent->print();
  TDSobj = TDSdir + "/Digi";
  // Register a Digi:
  log << MSG::DEBUG << "current digiEvent pointer: " << digiEvent << endreq;
  sc = m_dataSvc->registerObject(TDSobj, digiEvent);
  if ( sc.isFailure() ) {
    log << MSG::DEBUG << "failed to register " << TDSobj << endreq;
    return sc;
  }
  log << MSG::DEBUG <<" Digi event registered in "<<TDSobj<< endreq;
  return sc;
}


StatusCode AdfDigiAlg::finalize()
{
  return StatusCode::SUCCESS;
}
