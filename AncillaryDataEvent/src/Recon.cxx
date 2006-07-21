#include "AncillaryDataEvent/Recon.h"

using namespace AncillaryData;

Recon::Recon(AncillaryData::Digi *digiEvent)
{
  setEventNumber(digiEvent->getEventNumber());
  setSpillNumber(digiEvent->getSpillNumber());
  XCU=0;
  YCU=0;
  ZCU=0;
  E_rec=0; 
  E_corr=0;
  PhiIn=0;
  ThetaIn=0;
  Dphi=0;
  Egamma=0;
  EgammaErr=0;
}

void Recon::print()
{
  std::cout<< " Ancillary Recon Event: "<<getEventNumber()<<" Spill Number: "<<getSpillNumber()<<std::endl;
  std::cout<< " --- number of Tagger Clusters: "<<m_taggerClusterColl.size()<<std::endl;
  for(std::vector<TaggerCluster>::iterator pos=m_taggerClusterColl.begin(); pos!=m_taggerClusterColl.end(); ++pos)
    (*pos).print();
  std::cout<< " --- number of QDC Hits       : "<<m_qdcHitColl.size()<<std::endl;
  for(std::vector<QdcHit>::iterator        pos=       m_qdcHitColl.begin(); pos!=       m_qdcHitColl.end(); ++pos)
    (*pos).print();
  std::cout<< " --- number of Scaler Hits       : "<<m_scalerHitColl.size()<<std::endl;
  for(std::vector<ScalerHit>::iterator        pos=       m_scalerHitColl.begin(); pos!=       m_scalerHitColl.end(); ++pos)
    (*pos).print();

}

std::vector<TaggerCluster> Recon::GetHighestClusters()
{
  SortClusters();
  if(m_taggerClusterColl.size()<=1) 
    return m_taggerClusterColl;
  std::vector<TaggerCluster> HigestsClusters;
  std::vector<TaggerCluster>::iterator pos=m_taggerClusterColl.begin();
  TaggerCluster selectedCluster=(*pos);
  pos++;
  int i=1;
  while(pos<m_taggerClusterColl.end())
    {
      TaggerCluster newCluster=(*pos);
      if(selectedCluster.getModuleId()==newCluster.getModuleId() && 
	 selectedCluster.getLayerId()==newCluster.getLayerId())
	{
	  selectedCluster=MaxCluster(newCluster,selectedCluster);
	}
      else
	{
	  HigestsClusters.push_back(selectedCluster);
	  selectedCluster=newCluster;
	}
      pos++;
    }
  return HigestsClusters;
}

void Recon::SortClusters()
{
  std::sort(m_taggerClusterColl.begin(),m_taggerClusterColl.end(),ClusterSortPredicate);
}


void Recon::computePositions(AncillaryGeometry *geometry)
{
  std::vector<TaggerCluster> higestClusters = GetHighestClusters();
  for (unsigned int i = 0 ; i < N_MODULES ; i++)
    {
      X[i] = 0.0 ; Y[i] = 0.0; Z[i] = 0.0;
    }
  
  for(std::vector<TaggerCluster>::iterator pos=higestClusters.begin();pos!=higestClusters.end(); ++pos)
    {
      (*pos).computePosition(geometry);
      const unsigned int M= (*pos).getModuleId();
      const unsigned int L= (*pos).getLayerId();
      unsigned int View = geometry->getView(M); // 1 is Y (strip along Z) , 0 is Z (strip along Y)
      
      X[M] = geometry->getX(M);
      if(View==0)
	{
	  if(L==0)
	    {
	      Y[M] = (*pos).getY();
	    }
	  else 
	    {
	      Z[M] = (*pos).getZ();
	    }
	}
      else
	{
	  if(L==1)
	    {
	      Y[M] = (*pos).getY();
	    }
	  else 
	    {
	      Z[M] = (*pos).getZ();
	    }
	}
    }
}

void Recon::reconstructEnergy(AncillaryGeometry *geometry)
{
   E_rec  = 0.0;
   E_corr =0.0;
   double bt=geometry->getBL();
   /*
     double alpha = atan2(positiony[1]-positiony[0]),dist1);
     double beta  = atan2(positiony[3]-positiony[2]),dist2);
     if(fabs(sin(beta)-sin(alpha))>0)
     E_rec = 300.*bt/(sin(beta)-sin(alpha));
 vv0[0] = dist1;
    vv0[1] = positiony[1]-positiony[0];
    vv0[2] = positionz[1]-positionz[0];
    mynorm = sqrt(vv0[0]*vv0[0]+vv0[1]*vv0[1]+vv0[2]*vv0[2]);
    vv0[0] /= mynorm;
    vv0[1] /= mynorm;
    vv0[2] /= mynorm;
    vv1[0] = dist2;
    vv1[1] = positiony[3]-positiony[2];
    vv1[2] = positionz[3]-positionz[2];
    mynorm = sqrt(vv1[0]*vv1[0]+vv1[1]*vv1[1]+vv1[2]*vv1[2]);
    vv1[0] /= mynorm;
    vv1[1] /= mynorm;
    vv1[2] /= mynorm;
    alpha2 = ( asin(vv0[2]) + asin(vv1[2]))/2;
    dir_rec_x = vv0[0];
    dir_rec_y = vv0[1];
    dir_rec_z = vv0[2];
    if(fabs(cos(alpha2))>0)
    E_rec2 /= cos(alpha2);
   */
}  

/*
  AD_TIMESTAMP 	  Trigger Time-stamp measured in the AD system
  AD_PID 	 particle ID code from AD data (use same as G4 particle code?)
  TAG_PHI_IN 	 electron incoming angle before magnet measured by silicon chambers 1 and 2 in the bending plane
  TAG_THETA_IN 	 electron incoming angle before magnet measured by silicon chambers 1 and 2 in the bending plane
  TAG_XYZ[3,4] 	 X,Y,Z coordinates of highest cluster in each of the four silicon tagger station
  TAG_XYZ_IN_CU[3] 	 X,Y,Z coordinates of photon impact point on CU (assuming photon collinear with e beam)
  TAG_DPHI 	 electron deflection angle after magnet in the bending plane
  TAG_EGAMMA 	 photon energy (Ebeam - Edeflected_electron)
  TAG_EGAMMA_ERR 	 error on photon energy measurement
  CRNKV_PHA[2] 	 pulse height amplitude in the 2 cerenkov
  SCINT_PHA[8] 	 array of floats - PHA for scintillators in the setup (8 is a placeholder for a reasonable number, might change)
  TRD_NUM[16]

  void Recon::computeFinalTupla()
  {
  
  }
*/
