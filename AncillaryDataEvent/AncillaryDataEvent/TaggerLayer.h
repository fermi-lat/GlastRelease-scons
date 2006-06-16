#ifndef TAGGERLAYER_HH
#define TAGGERLAYER_HH

// Standard stuff.
#include <vector>

// ROOT stuff.
#include <TH1F.h>

// Our own classes.
#include "TaggerCluster.h"


class TaggerLayer

{
 public:
 
  TaggerLayer(int layerId, int moduleId);
  ~TaggerLayer();

  void bookHistograms();
  void processData();
  void reset();
  //void findHits();
  void findClusters();
  void fillHistograms();
  void writeHistogramsToFile();

  std::vector<TaggerHit*> m_hits;  
  int            getId()                        const {return m_layerId;}  
  TaggerChannel *getChannel(int channelId)      const {return m_channels[channelId];}
  TaggerHit     *getHit(int hitId)              const {return m_hits[hitId];}
  unsigned int   getNumberOfHits()              const {return m_hits.size();}
  unsigned int   getNumberOfClusters()          const {return m_clusters.size();}
  TH1F          *getNumberOfHitsHistogram()     const {return m_numberOfHitsHistogram;}
  TH1F          *getNumberOfClustersHistogram() const {return m_numberOfClustersHistogram;}
  TH1F          *getClusterSizeHistogram()      const {return m_clusterSizeHistogram;}
  TaggerCluster *getHighestCluster()            const {return m_highestCluster;}

 private:
  
  bool m_DEBUG;
  int  m_layerId;
  int  m_moduleId;
  TaggerChannel               *m_channels[N_CHANNELS_PER_LAYER];
  std::vector<TaggerCluster*>  m_clusters;
  TaggerCluster               *m_highestCluster;
  TH1F                        *m_numberOfHitsHistogram;
  TH1F                        *m_hitPositionHistogram;
  TH1F                        *m_hitPulseHeightHistogram;
  TH1F                        *m_numberOfClustersHistogram;
  TH1F                        *m_clusterSizeHistogram;
  TH1F                        *m_clusterPositionHistogram;
  TH1F                        *m_clusterPulseHeightHistogram;
  TH1F                        *m_highestClusterSizeHistogram;
  TH1F                        *m_highestClusterPositionHistogram;
  TH1F                        *m_highestClusterPulseHeightHistogram;
  TH1F                        *m_clusterSNRatioHistogram;
  TH1F                        *m_clusterHighestHitSNRatioHistogram;
  TH1F                        *m_highestClusterSNRatioHistogram;
  TH1F                        *m_highestClusterHighestHitSNRatioHistogram;
  TH1F                        *m_clusterEtaHistogram;
  TH1F                        *m_highestClusterEtaHistogram;
};
#endif
