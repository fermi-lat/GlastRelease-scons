#ifndef TAGGERCLUSTER_HH
#define TAGGERCLUSTER_HH

// Standard stuff.
#include <vector>
#include <math.h>

// Our own classes.
#include "TaggerHit.h"

class TaggerCluster

{
 public:
 
  TaggerCluster();
  ~TaggerCluster();
  void append(TaggerHit *hit);
  void calculateProperties();
  int  compareHits(TaggerHit *hit1, TaggerHit *hit2);
   
  unsigned int getSize()              const {return m_hits.size();}
  double       getPosition()          const {return m_baricenterPosition;}
  double       getPulseHeight()       const {return m_totalPulseHeight;}
  double       getNoise()             const {return m_totalNoise;}
  TaggerHit    *getHit(int hitId)     const {return m_hits[hitId];}
  std::vector<TaggerHit*> getHits()   const {return m_hits;}
  double       getSNRatio()           const {return m_signalToNoiseRatio;}
  double       getHighestHitSNRatio() const {return m_highestHitSignalToNoiseRatio;}
  double       getEta()               const {return m_eta;}
 
 private:
  
  bool       m_DEBUG;
  TaggerHit *m_highestHit;
  double     m_baricenterPosition;
  double     m_totalPulseHeight;
  double     m_totalNoise;
  double     m_signalToNoiseRatio;
  double     m_highestHitSignalToNoiseRatio;
  double     m_eta;
  std::vector<TaggerHit*> m_hits;
  
};
#endif
