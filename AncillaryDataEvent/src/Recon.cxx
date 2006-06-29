#include "AncillaryDataEvent/Recon.h"
using namespace AncillaryData;
Recon::Recon(AncillaryData::Digi *digiEvent)
{
  setEventNumber(digiEvent->getEventNumber());
  setSpillNumber(digiEvent->getSpillNumber());
}

void Recon::print()
{
  std::cout<< " Ancillary Recon Event: "<<getEventNumber()<<" Spill Number: "<<getSpillNumber()<<std::endl;
  std::cout<< " --- number of Tagger Clusters: "<<m_taggerClusterColl.size()<<std::endl;
  for(std::vector<TaggerCluster>::iterator pos=m_taggerClusterColl.begin();pos!=m_taggerClusterColl.end(); ++pos)
    (*pos).print();
  std::cout<< " --- number of QDC Hits       : "<<m_qdcHitColl.size()<<std::endl;
  for(std::vector<QdcHit>::iterator pos=m_qdcHitColl.begin();pos!=m_qdcHitColl.end(); ++pos)
    (*pos).print();

}

