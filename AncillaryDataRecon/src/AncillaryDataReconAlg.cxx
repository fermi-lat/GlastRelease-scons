// Gaudi specific include files
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "facilities/Util.h"

#include "AncillaryDataEvent/Digi.h"

//#include "AncillaryDataUtil/AncillaryDataServer.h"

const std::string TDSdir = "/Event/AncillaryEvent";

class AncillaryDataReconAlg : public Algorithm
{
  
public:
  AncillaryDataReconAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~AncillaryDataReconAlg() {}
  
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();
  
  
  AncillaryData::Digi *GetDigi();
  StatusCode RegisterRecon();
  StatusCode MakeClusters();
  StatusCode Recon();
private:
  IDataProviderSvc    *m_dataSvc;
};

static const AlgFactory<AncillaryDataReconAlg>  Factory;
const IAlgFactory& AncillaryDataReconAlgFactory = Factory;

AncillaryDataReconAlg::AncillaryDataReconAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : Algorithm(name, pSvcLocator) 
{
  // Input parameters that may be set via the jobOptions file
  //  declareProperty("dataFilePath",dataFilePath="$(ADFREADERROOT)/data/CR_DAQBARI_330000723.bin");
}


StatusCode AncillaryDataReconAlg::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << "Initialize AncillaryDataReconAlg"<< endreq;
  setProperties();
  IService* iService = 0;
  sc = serviceLocator()->service("EventDataSvc", iService, true);
  
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "could not find DataSvc !" << endreq;
    return sc;
  }
  m_dataSvc = dynamic_cast<IDataProviderSvc*>(iService);
  //DEFAULT ARGS:
  return sc;
}
StatusCode AncillaryDataReconAlg::execute()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  AncillaryData::Digi *digiEvent = GetDigi();
  log << MSG::DEBUG << "digiEvent "<< digiEvent << endreq;
  if(digiEvent) digiEvent->print();
  return sc;
}


StatusCode AncillaryDataReconAlg::finalize()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  return sc;
}

AncillaryData::Digi *AncillaryDataReconAlg::GetDigi()
{
  int I=0;
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name()); 
  // Get the ADF event:
  std::string TDSobj = TDSdir + "/Digi";
  log << MSG::DEBUG << " Get "<<TDSobj<<" from the TDS " << endreq;
  //////////////////////////////////////////////////
  // 1) test :create a sub dir (got from BariMcToolHitCol):
  log << MSG::DEBUG << " GetDigi: ...create " << TDSdir << endreq;
  DataObject* pNode = 0;
  sc = m_dataSvc->retrieveObject(TDSdir, pNode);
  if ( sc.isFailure() ) 
    {
      sc = m_dataSvc->registerObject(TDSdir, new DataObject);
      log << MSG::DEBUG << " RegisterDigi: registering " << TDSdir << endreq;
      if( sc.isFailure() ) 
	{
	  log << MSG::ERROR << "could not register " << TDSdir << endreq;
	  return 0;
	}
    }
  SmartDataPtr<AncillaryData::Digi> currentEvent(eventSvc(),TDSobj);
  
  if( !currentEvent) 
    {
      log<<MSG::WARNING<<"AncillaryDataReconAlg::GetDigi : empty TDS!" <<endreq;
      return 0;
    }
  log << MSG::DEBUG << "So far so good! " << endreq;
  return currentEvent;
  
}

StatusCode AncillaryDataReconAlg::RegisterRecon()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  return sc;
}

StatusCode AncillaryDataReconAlg:: MakeClusters()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  return sc;
}

StatusCode AncillaryDataReconAlg::Recon()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  return sc;
}

