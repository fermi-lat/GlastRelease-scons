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

  StatusCode RegisterTDSDir();
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
  declareProperty("dataFilePath",dataFilePath="$(ADFREADERDATAPATH)/CR_01_v3.bin");
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
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "initialize" << endreq;
  sc = RegisterTDSDir();
  sc = RegisterAdf();
  if (sc.isFailure()) {
    log << MSG::INFO << "Failed to register ADF" << endreq;
    return sc;
  }
  sc = RegisterDigi();
  if (sc.isFailure()) {
    log << MSG::INFO << "Failed to register Digi" << endreq;
    return sc;
  }
  return sc;
}

StatusCode AdfReaderAlg::RegisterTDSDir()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << "RegisterTDSDir " << endreq;
  // 1) REGISTER THE DIRECTORY:
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
    }
  return sc;
}

StatusCode AdfReaderAlg::RegisterAdf()
{  
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << "RegisterAdf " << endreq;
  //////////////////////////////////////////////////
  AncillaryData::AdfEvent *currentEvent = m_dataServer->getNextEvent();
  std::string TDSobj = TDSdir + "/AdfEvent";
  // 2) REGISTER THE OBJECT:
  sc = m_dataSvc->registerObject(TDSobj, currentEvent);
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "failed to register " << TDSobj << endreq;
    return sc;
  }
  log << MSG::DEBUG <<" Event registered in "<<TDSobj<< endreq;
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
  sc = m_dataSvc->registerObject(TDSobj, digiEvent);
  if ( sc.isFailure() ) {
    log << MSG::DEBUG << "failed to register " << TDSobj << endreq;
    return sc;
  }
  log << MSG::DEBUG <<" Digi event registered in "<<TDSobj<< endreq;
  return sc;
}


StatusCode AdfReaderAlg::finalize()
{
  std::cout<<" Finalize AdfReaderAlg "<<std::endl;
  m_dataServer->close();
  //  m_tagger->writeOutputRootFile(outputFilePath);

  return StatusCode::SUCCESS;
}
