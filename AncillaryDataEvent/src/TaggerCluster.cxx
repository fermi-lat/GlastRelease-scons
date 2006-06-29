#include "AncillaryDataEvent/TaggerCluster.h"
using namespace AncillaryData;

void TaggerCluster::calculateProperties()
{  
  m_properties=true;
  //sort(m_hits, m_hits.size(), sizeof(TaggerHit), compareHits);
  std::vector<TaggerHit>::iterator hitIterator;
  // Evaluate baricenter position, highest strip and total pulse height.
  double highestPulseHeight = 0.0;
  for (hitIterator = m_hits.begin(); hitIterator != m_hits.end(); hitIterator++){
    if ((*hitIterator).getPulseHeight() > highestPulseHeight){
      m_highestHit = (*hitIterator);
      highestPulseHeight = (*hitIterator).getPulseHeight();
    }
    m_baricenterPosition += (*hitIterator).getStripId()*(*hitIterator).getPulseHeight();
    m_totalPulseHeight   += (*hitIterator).getPulseHeight();
    //    m_totalNoise         += pow(((*hitIterator).getNoise()), 2.0);
  }
  m_baricenterPosition /= m_totalPulseHeight;
  m_baricenterPosition *= STRIPS_PITCH;
  //  m_totalNoise          = sqrt(m_totalNoise);
  // Evaluate signal to noise ratio.
  //  m_signalToNoiseRatio           = m_totalPulseHeight/m_totalNoise;
  //  m_highestHitSignalToNoiseRatio = m_highestHit.getPulseHeight()/m_highestHit.getNoise();
  // Evaluate Eta.
  double nextHighestPulseHeight = 0.0;
  if (getSize() > 1){
    for (hitIterator = m_hits.begin(); hitIterator != m_hits.end(); hitIterator++){
      if((abs((*hitIterator).getStripId() - m_highestHit.getStripId()) == 1) &&
         ((*hitIterator).getPulseHeight() > nextHighestPulseHeight)){
        nextHighestPulseHeight = (*hitIterator).getPulseHeight();
      }
    }
    m_eta = (highestPulseHeight - nextHighestPulseHeight)/(highestPulseHeight + nextHighestPulseHeight);
  }
}

void TaggerCluster::print()
{
  if(!m_properties) calculateProperties();
  std::cout<<" Size: "<<getSize()<<" Baricenter: "<<getPosition()<<" PH: "<<getPulseHeight()<<" Eta: "<<getEta()<<std::endl;
  for (std::vector<TaggerHit>::iterator hitIterator = m_hits.begin(); hitIterator != m_hits.end(); hitIterator++)
    (*hitIterator).print();
}
