#ifndef DataStructures_H
#define DataStructures_H


#include "Event/Recon/CalRecon/CalXtalRecData.h"

#include <vector>


/**
 * Class GcrHit
 * 
 */
class GcrHit:virtual public ContainedObject {

	public:
	
	  /**GcrHit(Event::CalXtalRecData* xtalData, bool isGoodHit, double correctedEnergy){
		  m_xtalData = xtalData;
		  m_isGoodHit = isGoodHit;
		  m_correctedEnergy = correctedEnergy;
	  }*/

	  GcrHit(Event::CalXtalRecData* xtalData, bool isGoodHit, double correctedEnergy){
		  m_xtalData = xtalData;
		  m_isGoodHit = isGoodHit;
		  m_correctedEnergy = correctedEnergy;
	  }


	  Event::CalXtalRecData* getXtalData(){return m_xtalData;}

	  void setXtalData(Event::CalXtalRecData* xtalData){m_xtalData = xtalData;}

	  bool getIsGoodHit(){return m_isGoodHit;}

	  void setIsGoodHit(bool isGoodHit){m_isGoodHit = isGoodHit;}

	  double getCorrectedEnergy(){return m_correctedEnergy;}

	  void setCorrectedEnergy(double correctedEnergy){m_correctedEnergy = correctedEnergy; }


	private:

	   Event::CalXtalRecData* m_xtalData;
	   bool m_isGoodHit;
	   double m_correctedEnergy;

};

// Define a vector of GcrHits
typedef std::vector<GcrHit> GcrHitsVec;

/**
 * Class GcrCluster
 * 
 */
class GcrCluster:virtual public ContainedObject {
public: 

   GcrCluster(const GcrHitsVec& hitsVec, double totalCorrectedEnergy, bool isGoodCluster){
	   m_hitsVec = hitsVec;
	   m_totalCorrectedEnergy = totalCorrectedEnergy;
	   m_isGoodCluster = isGoodCluster;
	   m_totalRawEnergy = -900;
	   m_totalPathLength = -900;
   }

   GcrCluster(double totalCorrectedEnergy, bool isGoodCluster){
	   m_totalCorrectedEnergy = totalCorrectedEnergy;
	   m_isGoodCluster = isGoodCluster;
	   m_totalRawEnergy = -900;
	   m_totalPathLength = -900;
   }

   GcrCluster(const GcrHitsVec& hitsVec, double totalCorrectedEnergy, bool isGoodCluster, double totalRawEnergy, double totalPathLength){
	   m_hitsVec = hitsVec;
	   m_totalCorrectedEnergy = totalCorrectedEnergy;
	   m_isGoodCluster = isGoodCluster;
	   m_totalRawEnergy = totalRawEnergy;
	   m_totalPathLength = totalPathLength;
   }

   GcrCluster(double totalCorrectedEnergy, bool isGoodCluster, double totalRawEnergy, double totalPathLength){
	   m_totalCorrectedEnergy = totalCorrectedEnergy;
	   m_isGoodCluster = isGoodCluster;
	   m_totalRawEnergy = totalRawEnergy;
	   m_totalPathLength = totalPathLength;
   }

    virtual ~GcrCluster() {
        m_hitsVec.clear();
    }

   double getTotalCorrectedEnergy(){return m_totalCorrectedEnergy;}
   void setTotalCorrectedEnergy(double totalCorrectedEnergy){m_totalCorrectedEnergy = totalCorrectedEnergy;}

   double getTotalRawEnergy(){return m_totalRawEnergy;}
   void setTotalRawEnergy(double totalRawEnergy){m_totalRawEnergy = totalRawEnergy;}

   double getTotalPathLength(){return m_totalPathLength;}
   void setTotalPathLength(double totalPathLength){m_totalPathLength = totalPathLength;}


   GcrHitsVec& getHitsVec(){return m_hitsVec;}
   void setHitsVec(const GcrHitsVec& hitsVec){m_hitsVec = hitsVec;}

   bool getIsGoodCluster(){return m_isGoodCluster;}
   void setIsGoodCluster(bool isGoodCluster){m_isGoodCluster = isGoodCluster;}

   int getMultiplicity(){int mult=0; return mult;}//method not implemented yet

private:

   GcrHitsVec m_hitsVec;
   float m_totalCorrectedEnergy;
   bool m_isGoodCluster;
   double m_totalRawEnergy;
   double m_totalPathLength;

};

// Define a vector of GcrClusters
typedef std::vector<GcrCluster> GcrClustersVec;

/**
 * Class GcrLayer
 * 
 */

class GcrLayer:virtual public ContainedObject{

public:

  GcrLayer(int layerNb, GcrClustersVec* gcrClustersVec, bool isGoodLayer, bool matchEneCrit){
	m_layerNb = layerNb;
	m_clustersVec = gcrClustersVec;
	m_isGoodLayer = isGoodLayer;
	m_matchEneCrit = matchEneCrit;
  }
  
  GcrClustersVec* getClustersVec(){return m_clustersVec;}

  void setClustersVec(GcrClustersVec* clustersVec){m_clustersVec = clustersVec;}

  bool getIsGoodLayer(){return m_isGoodLayer;}

  void setIsGoodLayer(bool value){m_isGoodLayer = value;}

  bool getMatchEneCrit(){return m_matchEneCrit;}

  void setMatchEneCrit(bool value){m_matchEneCrit = value;}

  int getLayerNb(){return m_layerNb;}
  
  void setLayerNb(int layerNb){m_layerNb = layerNb;}

private:

   int m_layerNb;
   GcrClustersVec* m_clustersVec;
   bool m_isGoodLayer;
   bool m_matchEneCrit;
   

};

// Define a vector of GcrLayers
typedef std::vector<GcrLayer> GcrLayersVec;


/**
 * Class GcrTower
 * 
 */
 
 class GcrTower:virtual public ContainedObject{

public:

  //GcrTower(int towerNb, GcrLayersVec* gcrLayersVec){
//	m_towerNb = towerNb;
//	m_layersVec = gcrLayersVec;
//  }

  GcrTower(int towerNb){
	m_towerNb = towerNb;
  }

  virtual ~GcrTower() {
      m_layersVec.clear();
  }

  GcrLayersVec& getLayersVec(){return m_layersVec;}
  void setLayersVec(const GcrLayersVec& layersVec){m_layersVec = layersVec;}

  int getTowerNb(){return m_towerNb;} 
  void setTowerNb(int towerNb){m_towerNb = towerNb;}



private:

   int m_towerNb;
   GcrLayersVec m_layersVec;
  

};

// Define a vector of GcrTowers
typedef std::vector<GcrTower> GcrTowersVec;

/**
 * Class LayerSpEnData
 * 
 */
 
 class LayerSpEnData:virtual public ContainedObject{

public:

  LayerSpEnData(double ePeak, double sigma){
	m_ePeak = ePeak;
	m_sigma = sigma;
  }
  
  double getEPeak(){return m_ePeak;} 
  void setEPeak(double ePeak){m_ePeak = ePeak;}

  double getSigma(){return m_sigma;} 
  void setSigma(double sigma){m_sigma = sigma;}


private:

   double m_ePeak;
   double m_sigma;
  

};

// Define a vector of LayerSpEnData
static const int NLAY = 8;
typedef std::vector<LayerSpEnData> SpEnergiesVec;


/**
 * Class SpEnergiesInterval
 * 
 */
 
 class SpEnergiesInterval:virtual public ContainedObject{

public:

  SpEnergiesInterval(){
	m_minE = -9999.99;
	m_maxE = -9999.99;
  }
  
  SpEnergiesInterval(double minE, double maxE){
	m_minE = minE;
	m_maxE = maxE;
  }
  
  double getMinE(){return m_minE;} 
  void setMinE(double minE){m_minE = minE;}

  double getMaxE(){return m_maxE;} 
  void setMaxE(double maxE){m_maxE = maxE;}


private:

   double m_minE;
   double m_maxE;
  

};

 
  
 /**
 * Class ZvecEntry: for each Z, A vector containing Spected Deposited Energy and Sigma for each Layer is given 
 * 
 */
 
 class ZVectEntry:virtual public ContainedObject{

public:

  ZVectEntry(int z, SpEnergiesVec* spEnergiesVec){
	m_z = z;
	m_spEnergiesVec = spEnergiesVec;
  }
  
  int getZ(){return m_z;} 
  void setZ(int z){m_z = z;}

  SpEnergiesVec* getSpEnergiesVec(){return m_spEnergiesVec;}
  void setSpEnergiesVec(SpEnergiesVec* spEnergiesVec){m_spEnergiesVec = spEnergiesVec;}




private:

   int m_z;
   SpEnergiesVec* m_spEnergiesVec;
  

};



// Define a vector of ZvectEntries
/**typedef std::vector<ZVectEntry> ZVect;*/
 


  
 /**
 * Class ZvecEntry2:  for each Z, an interval of spected deposited energy is given, only for first layer.
 * 
 */
 
 class ZVectEntry2:virtual public ContainedObject{

public:

  ZVectEntry2(int z, SpEnergiesInterval* spEnergiesInterval){
	m_z = z;
	m_spEnergiesInterval = spEnergiesInterval;
  }
  
  int getZ(){return m_z;} 
  void setZ(int z){m_z = z;}

  SpEnergiesInterval* getSpEnergiesInterval(){return m_spEnergiesInterval;}
  void setSpEnergiesInterval(SpEnergiesInterval* spEnergiesInterval){m_spEnergiesInterval = spEnergiesInterval;}




private:

   int m_z;
   SpEnergiesInterval* m_spEnergiesInterval;
  

};

typedef std::vector<ZVectEntry2> ZVect;





 

#endif //DataStructures_H
