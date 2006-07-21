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

#include "AncillaryDataUtil/AncillaryGeometry.h"

#include "CalibData/CalibModel.h"
#include "CalibData/Anc/AncCalibTaggerGain.h"
#include "CalibData/Anc/AncCalibTaggerPed.h"
#include "CalibData/Anc/AncCalibQdcPed.h"


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

  StatusCode SubtractPedestal(AncillaryData::Digi *digiEvent);
  StatusCode RegisterRecon(AncillaryData::Recon *reconEvent);
  StatusCode MakeClusters(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent);
  StatusCode taggerRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent);
  StatusCode QdcRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent);
  StatusCode ScalerRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent);
private:
  IDataProviderSvc    *m_dataSvc;
  IDataProviderSvc    *m_pCalibDataSvc;
  AncillaryData::AncillaryGeometry   *m_geometry;
  bool m_SubtractPedestals;
  std::string m_geometryFilePath;
};

static const AlgFactory<AncillaryDataReconAlg>  Factory;
const IAlgFactory& AncillaryDataReconAlgFactory = Factory;

AncillaryDataReconAlg::AncillaryDataReconAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : Algorithm(name, pSvcLocator) 
{
  // Input parameters that may be set via the jobOptions file
  declareProperty("pedestalSubtraction",m_SubtractPedestals=true);
  declareProperty("geometryFilePath",m_geometryFilePath="$(ANCILLARYDATAUTILROOT)/data/Geometry_v0.txt");
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
  
  sc = service("CalibDataSvc", m_pCalibDataSvc, true);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get IDataProviderSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;  
  }
  facilities::Util::expandEnvVar(&m_geometryFilePath);
  log << MSG::INFO << "loading geometry from " << m_geometryFilePath << endreq;
  m_geometry = new AncillaryData::AncillaryGeometry(m_geometryFilePath);
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

      AncillaryData::Recon *reconEvent = new AncillaryData::Recon(digiEvent);
      if(m_SubtractPedestals)
	{
	  sc = SubtractPedestal(digiEvent);
	  std::cout<<"------------ DIGI EVENT SUBRACTED PED !:-----------"<<std::endl;
	  digiEvent->print();
	  std::cout<<"---------------------------------------------------"<<std::endl;
	}
      sc = taggerRecon(digiEvent,reconEvent);
      sc = QdcRecon(digiEvent,reconEvent);      
      sc = ScalerRecon(digiEvent,reconEvent);      
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
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name()); 
  // Get the ADF event:
  const std::string TDSobj = TDSdir + "/Digi";
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
  
  const std::string TDSobj = TDSdir + "/Recon";
  sc = m_dataSvc->registerObject(TDSobj, reconEvent);
  if ( sc.isFailure() ) {
    log << MSG::DEBUG << "failed to register " << TDSobj << endreq;
    return sc;
  }
  log << MSG::DEBUG <<" Recon event registered in "<<TDSobj<< endreq;
  return sc;
} 

StatusCode AncillaryDataReconAlg::SubtractPedestal(AncillaryData::Digi *digiEvent)
{
  using CalibData::AncTaggerPed;
  using CalibData::AncQdcPed;
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  //////////////////////////////////////////////////
  // TAGGER:
  //////////////////////////////////////////////////
  std::string fullPathTaggerPed = "/Calib/ANC_TaggerPed/vanilla";
  DataObject *pObjectTagger;
  m_pCalibDataSvc->retrieveObject(fullPathTaggerPed, pObjectTagger);
  CalibData::AncCalibTaggerPed* pTaggerPeds = 0;
  pTaggerPeds = dynamic_cast<CalibData::AncCalibTaggerPed *> (pObjectTagger);
  if (!pTaggerPeds) {
    log << MSG::ERROR << "Dynamic cast to AncCalibTaggerPed failed" << endreq;
    return StatusCode::FAILURE;
  }
  
  std::vector<AncillaryData::TaggerHit> taggerHitCol=digiEvent->getTaggerHitCol();
  std::vector<AncillaryData::TaggerHit>::iterator taggerHitColI;
  
  for(taggerHitColI=taggerHitCol.begin(); taggerHitColI!=taggerHitCol.end();++taggerHitColI)
    {
      const unsigned int iMod  = (*taggerHitColI).getModuleId();
      const unsigned int iLay  = (*taggerHitColI).getLayerId();
      const unsigned int iChan = (*taggerHitColI).getStripId();
      CalibData::RangeBase* pTaggerPed = pTaggerPeds->getChan(iMod, iLay, iChan);
      AncTaggerPed* pT = dynamic_cast<AncTaggerPed * >(pTaggerPed);
      if (!pT) {
	log << MSG::ERROR << "Dynamic cast to AncTaggerPed failed M: " <<iMod<<" L: "<<iLay<<" C: "<<iChan<<endreq;
	return StatusCode::FAILURE;
      }
      const double pedestalValue  = static_cast<double>(pT->getVal());    
      const double PHPS           = (*taggerHitColI).getPulseHeight()-pedestalValue;
      const double RMS            = static_cast<double>(pT->getRNoise());
      
      if (PHPS<0)
	log << MSG::ERROR << " Negative pedestal subtracted Pulse haight "<<PHPS<<endreq;
      
      log << MSG::DEBUG << "ped = " << pedestalValue << endreq;
      log << MSG::DEBUG << "rNoise = " << RMS << endreq;
      log << MSG::DEBUG << "sNoise = " << pT->getSNoise() << endreq;
      log << MSG::DEBUG << "isBad = " <<  pT->getIsBad() << endreq;
      log << MSG::DEBUG << " PH ped subtract: "<<PHPS<< " ("<<PHPS/RMS<<") sigma"<<endreq;
      (*taggerHitColI).SubctractPedestal(pedestalValue, RMS);
      log << MSG::INFO << " Done pedestal subtraction for Tagger M: " <<iMod<<" L: "<<iLay
	  <<" C: "<<iChan<<" PH: ("<<PHPS<<")"<<endreq;
    }
  digiEvent->setTaggerHitCol(taggerHitCol);
  
  //////////////////////////////////////////////////
  // QDC
  std::string fullPathQdcPed = "/Calib/ANC_QdcPed/vanilla";
  DataObject *pObjectQdc;
  m_pCalibDataSvc->retrieveObject(fullPathQdcPed, pObjectQdc);
  CalibData::AncCalibQdcPed* pQdcPeds = 0;
  pQdcPeds = dynamic_cast<CalibData::AncCalibQdcPed *> (pObjectQdc);
  if (!pQdcPeds) {
    log << MSG::ERROR << "Dynamic cast to AncCalibQdcPed failed" << endreq;
    return StatusCode::FAILURE;
  }
  
  std::vector<AncillaryData::QdcHit> qdcHitCol=digiEvent->getQdcHitCol();
  std::vector<AncillaryData::QdcHit>::iterator qdcHitColI;
  for(qdcHitColI=qdcHitCol.begin(); qdcHitColI!=qdcHitCol.end();++qdcHitColI)
    {
      const unsigned int iChan = (*qdcHitColI).getQdcChannel();
      const unsigned int iMod  = (*qdcHitColI).getQdcModule();
      CalibData::RangeBase* pQdcPed = pQdcPeds->getQdcChan(iMod, iChan);
      
      AncQdcPed* pQ = dynamic_cast<AncQdcPed * >(pQdcPed);
      if (!pQ) {
	log << MSG::ERROR << "Dynamic cast to AncQdcPed failed M: " <<iMod<<" C: "<<iChan<<endreq;
	return StatusCode::FAILURE;
      }
      const double pedestalValue = static_cast<double> (pQ->getVal());    
      const double PHPS          = (*qdcHitColI).getPulseHeight()-pedestalValue;
      const double RMS           = static_cast<double> (pQ->getRms());
      log << MSG::DEBUG << "   ped = " << pedestalValue << endreq;
      log << MSG::DEBUG << "   rms = " << RMS << endreq;
      log << MSG::DEBUG << " isBad = " << pQ->getIsBad() << endreq;
      log << MSG::DEBUG << " PH ped subtract:"<<PHPS<<" ("<<PHPS/RMS<<") sigma"<<endreq;
      (*qdcHitColI).SubctractPedestal(pedestalValue, RMS);
      log << MSG::INFO << " Done pedestal subtraction for QDC M: " <<iMod<<" C: "<<iChan<<endreq;
    }
  digiEvent->setQdcHitCol(qdcHitCol);
  return sc;
}

StatusCode AncillaryDataReconAlg:: MakeClusters(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent)
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG <<" Making clusters "<<endreq; 
  // Get the tagger hits collection:
  std::vector<AncillaryData::TaggerHit> taggerHits = digiEvent->getTaggerHitCol();
  if (taggerHits.size()==0) return sc;
  std::vector<AncillaryData::TaggerCluster> taggerClusters;
  //  std::vector<AncillaryData::TaggerCluster>::iterator taggerClustersIterator;
  unsigned int i=0;
  AncillaryData::TaggerHit nextHit=taggerHits[i];
  nextHit.print();
  

  AncillaryData::TaggerCluster aCluster;
  aCluster.append(nextHit);
  //  taggerClusters.push_back(aCluster);

  AncillaryData::TaggerHit lastHit=nextHit;
  for(i=1; i<taggerHits.size();i++)
    {
      nextHit=taggerHits[i];
      nextHit.print();
      const signed int stripDist = 
	static_cast<signed int> (lastHit.getStripId())-
	static_cast<signed int> (nextHit.getStripId());
      if(nextHit.getModuleId()==lastHit.getModuleId() && 
	 nextHit.getLayerId()==lastHit.getLayerId() &&
	 abs(stripDist)<= static_cast<signed int> (AncillaryData::MAX_CLUSTER_GAP))
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
  taggerClusters.push_back(aCluster);
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

StatusCode AncillaryDataReconAlg::ScalerRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent)
{
  
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  reconEvent->setScalerHitColl(digiEvent->getScalerHitCol());
  return sc;
}

StatusCode AncillaryDataReconAlg::taggerRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent)
{
  MakeClusters(digiEvent,reconEvent);
  std::cout<<" Final Position of the higest cluster: "<<std::endl;
  std::cout<<"\t I \t X \t Y \t Z "<<std::endl; 
  reconEvent->computePositions(m_geometry);
  for (int i = 0;i<8;i++)
    {
      double x=reconEvent->getX(i);
      double y=reconEvent->getY(i);
      double z=reconEvent->getZ(i);
      std::cout<<"\t "<<i<<" \t "<<x<<" \t "<<y<<" \t "<<z<<std::endl; 
    }
}



