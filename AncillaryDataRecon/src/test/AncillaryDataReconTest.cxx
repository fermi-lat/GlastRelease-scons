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

#include <TTree.h>
#include <TFile.h>

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
  //  nextHit.print();
  
  AncillaryData::TaggerCluster aCluster;
  aCluster.append(nextHit);
  
  AncillaryData::TaggerHit lastHit=nextHit;
  for(i=1; i<taggerHits.size();i++)
    {
      nextHit=taggerHits[i];
      //      nextHit.print();
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

int main(int argc, char **argv)
{  
  std::string arg_name("");
  std::string dataFilePath="$(ADFREADERDATAPATH)/CR_01_v3.bin";
  int current_arg = 1;
  int Nevents=10;
  while(current_arg < argc)
    {
      arg_name = argv[current_arg];
      if( arg_name == "-n")
	{
	  Nevents=atoi(argv[++current_arg]);
	}
      else if( arg_name == "-f")
	{
	  dataFilePath=argv[++current_arg];
	}
      current_arg++;
    }
      
     
  
  // INITIALIZE THE GEOMETRY:
  std::string m_geometryFilePath="$(ANCILLARYDATAUTILDATAPATH)/Geometry_v0.dat";
  facilities::Util::expandEnvVar(&m_geometryFilePath); 
  std::cout<< "loading geometry from " << m_geometryFilePath << std::endl;
  AncillaryData::AncillaryGeometry   *m_geometry = new AncillaryData::AncillaryGeometry(m_geometryFilePath);
  //////////////////////////////////////////////////
  
  facilities::Util::expandEnvVar(&dataFilePath);
  AncillaryDataServer *m_dataServer = new AncillaryDataServer(dataFilePath);
  m_dataServer->open();
  m_dataServer->printInformation();
  //////////////////////////////////////////////////
  double X[4],Y[4],Z[4];
  double PX, PY, PZ;
  double E_rec, E_corr;
  double PhiIn,PhiOut;
  double Theta;
  unsigned int   evNumber,spillNumber;
  TTree* myEvent = new TTree("Event","Event");		
  myEvent->Branch("evNumber",&evNumber,"evNumber/i");
  myEvent->Branch("spillNumber",&spillNumber,"spillNumber/i");
  myEvent->Branch("X",&X,"X[4]/D");
  myEvent->Branch("Y",&Y,"Y[4]/D");
  myEvent->Branch("Z",&Z,"Z[4]/D");
  myEvent->Branch("PX",&PX,"PX/D");
  myEvent->Branch("PY",&PY,"PY/D");
  myEvent->Branch("PZ",&PZ,"PZ/D");
  myEvent->Branch("PhiIn",&PhiIn,"PhiIn/D");
  myEvent->Branch("PhiOut",&PhiOut,"PhiOut/D");
  myEvent->Branch("Theta",&Theta,"Theta/D");
  myEvent->Branch("E_rec",&E_rec,"E_rec/D");
  myEvent->Branch("E_corr",&E_corr,"E_corr/D");
  
  for (int i =0; i<Nevents; i++)
    {
      AncillaryData::AdfEvent *currentEvent = m_dataServer->getNextEvent();
      //      currentEvent->print();
      //      currentEvent->dump();
      //////////////////////////////////////////////////
      AncillaryData::Digi *digiEvent = new AncillaryData::Digi(currentEvent);
      digiEvent->print();
      //////////////////////////////////////////////////
      
      AncillaryData::Recon *reconEvent = new AncillaryData::Recon(digiEvent);
      // This is  in AncillaryDataReconAlg::taggerRecon(AncillaryData::Digi *digiEvent, AncillaryData::Recon *reconEvent)
      MakeClusters(digiEvent,reconEvent);
      reconEvent->ReconstructTagger(m_geometry);
      reconEvent->report();
      if(reconEvent->getNumberOfHigestClusters()==8)
	{
	  evNumber = reconEvent->getEventNumber();
	  spillNumber = reconEvent->getSpillNumber();
	  for (int i=0;i<4;i++)
	    {
	      X[i]    = reconEvent->getX(i);
	      Y[i]    = reconEvent->getY(i);
	      Z[i]    = reconEvent->getZ(i);
	    }
	  PX      = reconEvent->getPX();
	  PY      = reconEvent->getPY();
	  PZ      = reconEvent->getPZ();
	  PhiIn   = reconEvent->getPhiIn();
	  PhiOut  = reconEvent->getPhiOut();
	  Theta   = reconEvent->getTheta();
	  E_rec    = reconEvent->getReconstructedEnergy();
	  E_corr   = reconEvent->getCorrectedEnergy();
	  myEvent->Fill();
	}
    }
  TFile *myFile = new TFile("TestRecon.root","RECREATE");
  myEvent->Write();
  myFile->Close();
  m_dataServer->close();
  return 0;
}
