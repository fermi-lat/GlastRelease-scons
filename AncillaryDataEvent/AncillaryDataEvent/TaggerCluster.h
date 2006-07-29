#ifndef TAGGERCLUSTER_HH
#define TAGGERCLUSTER_HH

// Standard stuff.
#include <vector>
#include <math.h>

// Our own classes.
#include "TaggerHit.h"
#include "AdfEvent/TaggerParameters.h"
#include "AncillaryDataUtil/AncillaryGeometry.h"

namespace AncillaryData
{

  class TaggerCluster  
    {
    public:
      
      TaggerCluster()
	: m_baricenterPosition(0),m_totalPulseHeight(0),m_eta(0),m_properties(false)
	{ }
      
      ~TaggerCluster(){erase();}

      void initialize(double pos, double pulseHeight, double eta, const TaggerHit &highHit) {
          m_baricenterPosition = pos;
          m_totalPulseHeight = pulseHeight;
          m_eta = eta;
          m_highestHit = highHit; 
          m_properties = true;
      }
      
      void append(TaggerHit hit) 
	{
	  m_properties=false;
	  m_hits.push_back(hit);
	}
      
      void erase()
	{
	  m_properties=false;
	  m_hits.erase(m_hits.begin(),m_hits.end());
	}
      void calculateProperties();
      //int  compareHits(TaggerHit *hit1, TaggerHit *hit2);
      
      unsigned int getSize()              const {return m_hits.size();}
      double       getPosition()          const {return m_baricenterPosition;}
      double       getPulseHeight()       const {return m_totalPulseHeight;}
      unsigned int getModuleId()          const {return getHighestHit().getModuleId();}
      unsigned int getLayerId()           const {return getHighestHit().getLayerId();}
      TaggerHit    getHighestHit()        {if(!m_properties) calculateProperties(); return m_highestHit;}
      TaggerHit    getHighestHit()        const { return m_highestHit;}
      //      double       getNoise()             const {return m_totalNoise;}
      TaggerHit    getHit(int hitId)      const {return m_hits[hitId];}
      std::vector<TaggerHit> getHits()    const {return m_hits;}
      //      double       getSNRatio()           const {return m_signalToNoiseRatio;}
      //      double       getHighestHitSNRatio() const {return m_highestHitSignalToNoiseRatio;}
      void computePosition(AncillaryGeometry *geometry);
      double       getEta()               
	const {
	if(!m_properties) 
	  calculateProperties(); 
	return m_eta;
      }
      void print();
      
      double getX(){return m_X;}
      double getY(){return m_Y;}
      double getZ(){return m_Z;}
      
    private:
      TaggerHit  m_highestHit;
      double     m_baricenterPosition;
      double     m_totalPulseHeight;
      //      double     m_totalNoise;
      //      double     m_signalToNoiseRatio;
      //      double     m_highestHitSignalToNoiseRatio;
      double     m_eta;
      bool m_properties;
      double m_X, m_Y, m_Z;
      
      std::vector<TaggerHit> m_hits;    
    };
  
  inline bool ClusterSortPredicate(const TaggerCluster& cl1, const TaggerCluster& cl2)
    {
      const unsigned int M1 = cl1.getModuleId();
      const unsigned int L1 = cl1.getLayerId();
      const unsigned int M2 = cl2.getModuleId();
      const unsigned int L2 = cl2.getLayerId();
      const unsigned int I1 = L1+2*M1;
      const unsigned int I2 = L2+2*M2;
      return  (I1<I2);
    }
  
  inline  TaggerCluster MaxCluster(const TaggerCluster& cl1, const TaggerCluster& cl2)
    {
      if (cl1.getPulseHeight() > cl2.getPulseHeight()) 
	return cl1;
      return cl2;
    }

}
#endif

