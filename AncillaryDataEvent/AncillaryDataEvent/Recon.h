#ifndef ANCILLARYDATAEVENT_RECON_H
#define ANCILLARYDATAEVENT_RECON_H

#include "TaggerCluster.h"
#include "TaggerHit.h"
#include "QdcHit.h"
#include "Digi.h"
#include <vector>
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"


static const CLID& CLID_AncillaryDataReconEvent = InterfaceID("AncillaryDataReconEvent", 1, 0);

namespace AncillaryData
{
  class Recon : public DataObject 
    {
    public:
      Recon(){;}
      Recon(AncillaryData::Digi *digiEvent);
      static const CLID& classID()       { return CLID_AncillaryDataReconEvent; }
      
      void setTaggerClusters(std::vector<TaggerCluster> taggerClusterColl) {m_taggerClusterColl=taggerClusterColl;}
      void setQdcHitColl(std::vector<QdcHit> qdcHitColl) {m_qdcHitColl = qdcHitColl;}
      
      std::vector<TaggerCluster> getTaggerClusters() { return m_taggerClusterColl;}
      void print();
      void setEventNumber(unsigned eventNumber) {m_eventNumber=eventNumber;}
      unsigned getEventNumber(){return m_eventNumber;}
      void setSpillNumber(unsigned spillNumber) {m_spillNumber=spillNumber;}
      unsigned getSpillNumber(){return m_spillNumber;}
    private:
      unsigned m_eventNumber;
      unsigned m_spillNumber;
      std::vector<TaggerCluster> m_taggerClusterColl;
      std::vector<QdcHit> m_qdcHitColl;
    }; 
}
#endif
