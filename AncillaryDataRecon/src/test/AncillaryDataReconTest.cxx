#include <iostream>

#include "facilities/Util.h"

#include "AdfEvent/AdfEvent.h"
#include "AncillaryDataUtil/AncillaryDataServer.h"

#include "AncillaryDataEvent/Digi.h"
#include "AncillaryDataEvent/Recon.h"
#include "AncillaryDataEvent/TaggerHit.h"
#include "AncillaryDataEvent/TaggerCluster.h"
#include "AncillaryDataEvent/QdcHit.h"
#include "AncillaryDataUtil/AncillaryGeometry.h"


//////////////////////////////////////////////////
using namespace  AncillaryData;
//////////////////////////////////////////////////

// do something useful another day
//using namespace AncillaryData
int MakeClusters(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent)
{
  StatusCode sc = StatusCode::SUCCESS;
  std::cout<<" Making clusters "<<std::endl;
  // Get the tagger hits collection:
  std::vector<AncillaryData::TaggerHit> taggerHits = digiEvent->getTaggerHitCol();
  if (taggerHits.size()==0) return sc;
  std::vector<AncillaryData::TaggerCluster> taggerClusters;
  unsigned int i=0;
  AncillaryData::TaggerHit nextHit=taggerHits[i];
  nextHit.print();
  
  AncillaryData::TaggerCluster aCluster;
  aCluster.append(nextHit);
  
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
  return 0;
}

int main()
{  
  // INITIALIZE THE GEOMETRY:
  std::string m_geometryFilePath="$(ANCILLARYDATAUTILROOT)/data/Geometry_v0.dat";
  facilities::Util::expandEnvVar(&m_geometryFilePath); 
  std::cout<< "loading geometry from " << m_geometryFilePath << std::endl;
  AncillaryData::AncillaryGeometry   *m_geometry = new AncillaryData::AncillaryGeometry(m_geometryFilePath);
  //////////////////////////////////////////////////
  std::string dataFilePath="$(ADFREADERROOT)/data/CR_01_v3.bin";
  facilities::Util::expandEnvVar(&dataFilePath);
  AncillaryDataServer *m_dataServer = new AncillaryDataServer(dataFilePath);
  m_dataServer->open();
  m_dataServer->printInformation();
  for (int i =0; i<100; i++)
    {
      AncillaryData::AdfEvent *currentEvent = m_dataServer->getNextEvent();
      currentEvent->print();
      //      currentEvent->dump();
      //////////////////////////////////////////////////
      AncillaryData::Digi *digiEvent = new AncillaryData::Digi(currentEvent);
      digiEvent->print();
      //////////////////////////////////////////////////
      
      AncillaryData::Recon *reconEvent = new AncillaryData::Recon(digiEvent);
      // This is  in AncillaryDataReconAlg::taggerRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent)
      MakeClusters(digiEvent,reconEvent);
      
      std::cout<<" Final Position of the higest cluster: "<<std::endl;
      std::cout<<" \t I \t X \t Y \t Z "<<std::endl; 
      reconEvent->ReconstructTagger(m_geometry);
      for (int i = 0;i<4;i++)
	{
	  double x=reconEvent->getX(i);
	  double y=reconEvent->getY(i);
	  double z=reconEvent->getZ(i);
	  std::cout<<i<<" \t "<<x<<" \t "<<y<<" \t "<<z<<std::endl; 
	}
      std::cout<<" Reconstructed Energy: "<<reconEvent->getReconstructedEnergy()<<" "<<reconEvent->getCorrectedEnergy()<<" "<<std::endl;
      // This is QDC recon:
      reconEvent->setQdcHitColl(digiEvent->getQdcHitCol());      
      reconEvent->setScalerHitColl(digiEvent->getScalerHitCol());      
      std::cout<<"########  RECON EVENT PRINT:   ##################"<<std::endl;
      reconEvent->print();
      std::cout<<"##################################################"<<std::endl;
    }
  m_dataServer->close();
  return 0;
}
