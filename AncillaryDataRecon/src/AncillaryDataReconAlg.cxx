// Gaudi specific include files
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "facilities/Util.h"

#include "AncillaryDataEvent/Digi.h"
#include "AncillaryDataEvent/Recon.h"
#include "AncillaryDataEvent/TaggerHit.h"
#include "AncillaryDataEvent/TaggerCluster.h"
#include "AncillaryDataEvent/QdcHit.h"

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
  StatusCode RegisterTDSDir();

  StatusCode RegisterRecon(AncillaryData::Recon *reconEvent);
  StatusCode MakeClusters(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent);
  StatusCode QdcRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent);
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
  if(digiEvent) 
    {
      digiEvent->print();
      AncillaryData::Recon *reconEvent = new AncillaryData::Recon(digiEvent);
      sc = MakeClusters(digiEvent,reconEvent);
      sc = QdcRecon(digiEvent,reconEvent);
      sc = RegisterRecon(reconEvent);
      reconEvent->print();
    }
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
  return currentEvent;
  
}

StatusCode AncillaryDataReconAlg::RegisterTDSDir()
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

StatusCode AncillaryDataReconAlg::RegisterRecon(AncillaryData::Recon *reconEvent)
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  
  std::string TDSobj = TDSdir + "/Recon";
  sc = m_dataSvc->registerObject(TDSobj, reconEvent);
  if ( sc.isFailure() ) {
    log << MSG::DEBUG << "failed to register " << TDSobj << endreq;
    return sc;
  }
  log << MSG::DEBUG <<" Recon event registered in "<<TDSobj<< endreq;
  return sc;
} 

StatusCode AncillaryDataReconAlg:: MakeClusters(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent)
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  // Get the tagger hits collection:
  std::vector<AncillaryData::TaggerHit> taggerHits = digiEvent->getTaggerHitCol();
  if (taggerHits.size()==0) return sc;
  std::vector<AncillaryData::TaggerCluster> taggerClusters;
  std::vector<AncillaryData::TaggerCluster>::iterator taggerClustersIterator;
  //  std::vector<AncillaryData::TaggerHit>::iterator pos=taggerHits.begin();  
  int i=0;
  AncillaryData::TaggerHit nextHit=taggerHits[i];
  nextHit.print();


  AncillaryData::TaggerCluster aCluster;
  aCluster.append(nextHit);
  taggerClusters.push_back(aCluster);

  AncillaryData::TaggerHit lastHit=nextHit;
  for(i=1; i<taggerHits.size();i++)
    {
      nextHit=taggerHits[i];
      nextHit.print();
      
      if(nextHit.getModuleId()==lastHit.getModuleId() && 
	 nextHit.getLayerId()==lastHit.getLayerId() &&
	 nextHit.getStripId()==lastHit.getStripId()+1)
	{
	  aCluster.append(nextHit);
	  lastHit=nextHit;
	}
      else 
	{
	  taggerClusters.push_back(aCluster);
	  aCluster.erase();
	  aCluster.append(nextHit);
	  lastHit=nextHit;
	}
    }  
  reconEvent->setTaggerClusters(taggerClusters);
  return sc;
}




StatusCode AncillaryDataReconAlg::QdcRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent)
{
  
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  reconEvent->setQdcHitColl(digiEvent->getQdcHitCol());
  return sc;
}

