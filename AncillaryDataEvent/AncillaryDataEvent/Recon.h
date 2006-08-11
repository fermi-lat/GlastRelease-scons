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
      
      void setTaggerClusters(std::vector<TaggerCluster> taggerClusterCol) 
	{
	  m_taggerClusterCol=taggerClusterCol;
	}
      void eraseTaggerClusters(){m_taggerClusterCol.erase(m_taggerClusterCol.begin(),m_taggerClusterCol.end());}
      void appendTaggerCluster(const TaggerCluster clus) { m_taggerClusterCol.push_back(clus); }

      void setQdcHitCol(std::vector<QdcHit> qdcHitCol) {m_qdcHitCol = qdcHitCol;}
      void appendQdcHit(const QdcHit qdc) { m_qdcHitCol.push_back(qdc); }
      
      void setScalerHitCol(std::vector<ScalerHit> scalerHitCol) {m_scalerHitCol = scalerHitCol;}
      void appendScalerHit(const ScalerHit scaler) { m_scalerHitCol.push_back(scaler); }
      
      void SortClusters();
      std::vector<TaggerCluster> GetHighestClusters();
      void ComputeClustersProperties();
      void computePositions(AncillaryData::AncillaryGeometry* geometry);
      void reconstructEnergy(AncillaryData::AncillaryGeometry *geometry);
      //////////////////////////////////////////////////
      void ReconstructTagger(AncillaryData::AncillaryGeometry *geometry)
	{
	  ComputeClustersProperties();
	  computePositions(geometry);
	  reconstructEnergy(geometry);
	}
      
	
      const std::vector<TaggerCluster> getTaggerClusters() const { return m_taggerClusterCol; }
      const std::vector<QdcHit> getQdcHitCol() const { return m_qdcHitCol; }
      const std::vector<ScalerHit> getScalerHitCol() const { return m_scalerHitCol; }

      void print();
      void report();
      
      void setEventNumber(unsigned eventNumber) {m_eventNumber=eventNumber;}
      unsigned int getEventNumber() const {return m_eventNumber;}
      void setSpillNumber(unsigned spillNumber) {m_spillNumber=spillNumber;}
      unsigned int getSpillNumber() const {return m_spillNumber;}      

      void setNumberOfHighestClusters(unsigned int n) { m_NumberHigestClusters = n; }
      unsigned int getNumberOfHigestClusters() const { return  m_NumberHigestClusters;}
      
      void setNumberOfClusters(unsigned int l, unsigned int m, unsigned int n) { m_NumberClusters[l][m] = n; }
      unsigned int getNumberOfClusters(unsigned int l, unsigned int m) const { return  m_NumberClusters[l][m];}
      
      void setPos(const double *xArr, const double *yArr, const double *zArr) {
	unsigned int i;
	for (i=0; i<N_MODULES; i++) {
	  X[i] = xArr[i];
	  Y[i] = yArr[i];
	  Z[i] = zArr[i];
	}
      }
      double getX(unsigned int Module) const {return X[Module];}
      double getY(unsigned int Module) const {return Y[Module];}
      double getZ(unsigned int Module) const {return Z[Module];}
      
      void setIntersection(double x, double y, double z) {
          PX = x; PY = y; PZ = z; }
      double getPX() const {return PX;}
      double getPY() const {return PY;}
      double getPZ() const {return PZ;}
      void setPhi(double in, double out) {
          PhiIn = in; PhiOut = out; }
      double getPhiIn() const {return PhiIn;}
      double getPhiOut() const {return PhiOut;}
      void setThetaPhi(double t, double p) {
          Theta = t; Dphi = p; }
      double getTheta() const {return Theta;}
      double getDeltaPhi() const {return Dphi;}
      
      void setEnergy(double erec, double ecorr) {
          E_rec = erec;  E_corr = ecorr; }
      double getReconstructedEnergy() const{ return E_rec;}
      double getCorrectedEnergy() const {return E_corr;}
      
      void setEnergyError(double err1, double err2)
	{
	  Error_E_rec  = err1;
	  Error_E_corr = err2;
	}
      double getErrorErec() const {return Error_E_rec;}
      double getErrorEcorr() const {return Error_E_corr;}

    private:
      unsigned m_eventNumber;
      unsigned m_spillNumber;
      std::vector<TaggerCluster> m_taggerClusterCol;
      std::vector<QdcHit> m_qdcHitCol;
      std::vector<ScalerHit> m_scalerHitCol;
      double X[N_MODULES],Y[N_MODULES],Z[N_MODULES];

      /// coordinates of the intersection point between the incoming direction 
      /// and the outgoing direction
      double PX, PY, PZ;

      double E_rec, E_corr;
      double Error_E_rec, Error_E_corr;
      double PhiIn,PhiOut;
      double Theta;
      double Dphi;
      unsigned int m_NumberHigestClusters;
      unsigned int m_NumberTotalClusters;
      unsigned int m_NumberClusters[2][N_MODULES];
    }; 
}
#endif
