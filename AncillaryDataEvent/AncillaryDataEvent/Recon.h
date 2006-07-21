#ifndef ANCILLARYDATAEVENT_RECON_H
#define ANCILLARYDATAEVENT_RECON_H

#include "TaggerCluster.h"
#include "TaggerHit.h"
#include "QdcHit.h"
#include "ScalerHit.h"

#include "Digi.h"
#include <vector>
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"
#include "AncillaryDataUtil/AncillaryGeometry.h"


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
      void appendTaggerCluster(const TaggerCluster clus) { m_taggerClusterColl.push_back(clus); }

      void setQdcHitColl(std::vector<QdcHit> qdcHitColl) {m_qdcHitColl = qdcHitColl;}
      void appendQdcHit(const QdcHit qdc) { m_qdcHitColl.push_back(qdc); }
      
      void setScalerHitColl(std::vector<ScalerHit> scalerHitColl) {m_scalerHitColl = scalerHitColl;}
      void appendScalerHit(const ScalerHit scaler) { m_scalerHitColl.push_back(scaler); }
      
      void SortClusters();
      std::vector<TaggerCluster> GetHighestClusters();
      void computePositions(AncillaryData::AncillaryGeometry* geometry);
      void reconstructEnergy(AncillaryData::AncillaryGeometry *geometry);
      //////////////////////////////////////////////////
      void ReconstructTagger(AncillaryData::AncillaryGeometry *geometry)
	{
	  computePositions(geometry);
	  reconstructEnergy(geometry);
	}
      
	
      const std::vector<TaggerCluster> getTaggerClusters() const { return m_taggerClusterColl; }

      const std::vector<QdcHit>& getQdcHitCol() const { return m_qdcHitColl; }

      void print();
      void setEventNumber(unsigned eventNumber) {m_eventNumber=eventNumber;}
      unsigned getEventNumber() const {return m_eventNumber;}
      void setSpillNumber(unsigned spillNumber) {m_spillNumber=spillNumber;}
      unsigned getSpillNumber() const {return m_spillNumber;}

      double getX(unsigned int Module) {return X[Module];}
      double getY(unsigned int Module) {return Y[Module];}
      double getZ(unsigned int Module) {return Z[Module];}
      
      double getXCU() {return XCU;}
      double getYCY() {return YCU;}
      double getZCU() {return ZCU;}
      double getPhiIn(){return PhiIn;}
      double getThetaIn(){return ThetaIn;}
      double getDeltaPhi(){return Dphi;}
      double getEgamma(){return Egamma;}
      double getEgammaErr(){return EgammaErr;}
      
      double getReconstructedEnergy(){return E_rec;}
      double getCorrectedEnergy(){return E_corr;}

    private:
      unsigned m_eventNumber;
      unsigned m_spillNumber;
      std::vector<TaggerCluster> m_taggerClusterColl;
      std::vector<QdcHit> m_qdcHitColl;
      std::vector<ScalerHit> m_scalerHitColl;
      
      double X[N_MODULES],Y[N_MODULES],Z[N_MODULES];
      double XCU,YCU,ZCU;
      double E_rec, E_corr;
      double PhiIn,ThetaIn;
      double Dphi,Egamma,EgammaErr;
    }; 
}
#endif
