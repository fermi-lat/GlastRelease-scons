// Gaudi specific include files
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "facilities/Util.h"

#include "AdfEvent/AdfEvent.h"
#include "AncillaryDataEvent/Digi.h"
#include "AncillaryDataUtil/AncillaryDataServer.h"

const std::string TDSdir = "/Event/AncillaryEvent";

class AdfReaderAlg : public Algorithm
{

public:
  AdfReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~AdfReaderAlg() {}

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  StatusCode RegisterDigi();
  StatusCode RegisterAdf();

private:
  IDataProviderSvc    *m_dataSvc;
  
  AncillaryDataServer *m_dataServer;
  

  std::string dataFilePath;
};

static const AlgFactory<AdfReaderAlg>  Factory;
const IAlgFactory& AdfReaderAlgFactory = Factory;

AdfReaderAlg::AdfReaderAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : Algorithm(name, pSvcLocator) 
{
  // Input parameters that may be set via the jobOptions file
  declareProperty("dataFilePath",dataFilePath="$(ADFREADERROOT)/data/CR_DAQBARI_330000723.bin");
}

StatusCode AdfReaderAlg::initialize()
{
 // Purpose and Method:  Called once before the run begins.
  // This method:
  // 1) locates the Event Data Server.
  // 2) creates the AncillaryDataServer
  
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "initialize" << endreq;
  setProperties();
  IService* iService = 0;
  sc = serviceLocator()->service("EventDataSvc", iService, true);
  
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "could not find DataSvc !" << endreq;
    return sc;
  }
  m_dataSvc = dynamic_cast<IDataProviderSvc*>(iService);
  //DEFAULT ARGS:
  facilities::Util::expandEnvVar(&dataFilePath);
  log << MSG::INFO << "loading events from " << dataFilePath << endreq;

  m_dataServer = new AncillaryDataServer(dataFilePath);
  m_dataServer->open();
  m_dataServer->printInformation();
  return sc;  
}


StatusCode AdfReaderAlg::execute()
{  
  RegisterAdf();
  RegisterDigi();
}

StatusCode AdfReaderAlg::RegisterAdf()
{  
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << "execute AdfReaderAlg " << endreq;
  
  AncillaryData::AdfEvent *currentEvent = m_dataServer->getNextEvent();
  
  if ( !currentEvent ) 
    {
      log << MSG::DEBUG << "current event NULL, exiting!" << endreq;
      return sc;
    }
  
  currentEvent->print();
  
  
  // CREATE A DIGI EVENT
  AncillaryData::Digi *digiEvent = new AncillaryData::Digi(currentEvent);
  digiEvent->print();
  // 1) test :create a sub dir (got from BariMcToolHitCol):
  log << MSG::DEBUG << " ...create " << TDSdir << endreq;
  DataObject* pNode = 0;
  sc = m_dataSvc->retrieveObject(TDSdir, pNode);
  
  if ( sc.isFailure() ) 
    {
      sc = m_dataSvc->registerObject(TDSdir, new DataObject);
      log << MSG::DEBUG << " registering " << TDSdir << endreq;
      if( sc.isFailure() ) 
	{
	  log << MSG::ERROR << "could not register " << TDSdir << endreq;
	  return sc;
	}
      else
	{
	  log << MSG::DEBUG << TDSdir << " registered " << endreq;
	}
    }
  const std::string TDSobj = TDSdir + "/AdfEvent";
  
  // 2) Check if a AdfReader is already stored:
  /*
    log << MSG::DEBUG << "Check if " << TDSobj << " is already stored" << endreq;
    DataObject* dobj=0;
    sc = m_dataSvc->retrieveObject(TDSobj, dobj);
    if( sc.isFailure() )
    log << MSG::DEBUG << "could not retrieve " << TDSobj << endreq;
    else
    log << MSG::DEBUG << TDSobj << " found!" << endreq;
  */
  // 3) Register an AdfEvent:
  log << MSG::DEBUG << "currentEvent pointer: " << currentEvent << endreq;
  //  return StatusCode::SUCCESS;
  sc = m_dataSvc->registerObject(TDSobj, currentEvent);
  if ( sc.isFailure() ) {
    log << MSG::DEBUG << "failed to register " << TDSobj << endreq;
    return sc;
  }
  //  log << MSG::DEBUG << "Ancillary Event: "<<currentEvent->getEventNumber()<<" processed and stored " << TDSobj << " on TDS" << endreq;
  
  return sc;
}

StatusCode AdfReaderAlg::RegisterDigi()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  // Get the ADF event:
  std::string TDSobj = TDSdir + "/AdfEvent";
  log << MSG::DEBUG << " Get "<<TDSobj<<" from the TDS " << endreq;
  //////////////////////////////////////////////////
  // 1) test :create a sub dir (got from BariMcToolHitCol):
  log << MSG::DEBUG << " RegisterDigi: ...create " << TDSdir << endreq;
  DataObject* pNode = 0;
  sc = m_dataSvc->retrieveObject(TDSdir, pNode);
  if ( sc.isFailure() ) 
    {
      sc = m_dataSvc->registerObject(TDSdir, new DataObject);
      log << MSG::DEBUG << " RegisterDigi: registering " << TDSdir << endreq;
      if( sc.isFailure() ) 
	{
	  log << MSG::ERROR << "could not register " << TDSdir << endreq;
	  return sc;
	}
      else
	{
	  log << MSG::DEBUG << TDSdir << "RegisterDigi: registered " << endreq;
	}
    }  
  //////////////////////////////////////////////////
  // 1) method:
  /*DataObject* pNode = 0;
    sc = m_dataSvc->retrieveObject(TDSobj, pNode);
    AncillaryData::AdfEvent *currentEvent = dynamic_cast<AncillaryData::AdfEvent*>(pNode);
  */
  // 2) Method:
  SmartDataPtr<AncillaryData::AdfEvent> currentEvent(eventSvc(),TDSobj);
  if( !currentEvent) 
    {
      log<<MSG::WARNING<<"AdfReaderAlg::RegisterDigi : empty TDS!" <<endreq;
      return sc;
    }
  
  // CREATE A DIGI EVENT
  AncillaryData::Digi *digiEvent = new AncillaryData::Digi(currentEvent);
  digiEvent->print();
  TDSobj = TDSdir + "/Digi";
  // Register a Digi:
  log << MSG::DEBUG << "current digiEvent pointer: " << digiEvent << endreq;
  sc = m_dataSvc->registerObject(TDSobj, currentEvent);
  if ( sc.isFailure() ) {
    log << MSG::DEBUG << "failed to register " << TDSobj << endreq;
    return sc;
  }
  return sc;
}


StatusCode AdfReaderAlg::finalize()
{
  std::cout<<" Finalize AdfReaderAlg "<<std::endl;
  m_dataServer->close();
  //  m_tagger->writeOutputRootFile(outputFilePath);

  return StatusCode::SUCCESS;
}
