#include "AncillaryDataEvent/Recon.h"
#include <algorithm>

//using namespace AncillaryData;

namespace AncillaryData {

Recon::Recon(AncillaryData::Digi *digiEvent)
{
  setEventNumber(digiEvent->getEventNumber());
  setSpillNumber(digiEvent->getSpillNumber());
  setQdcHitColl(digiEvent->getQdcHitCol());      
  setScalerHitColl(digiEvent->getScalerHitCol());      

  PX =0.0;  PY=0.0;  PZ=0.0;
  E_rec=0; 
  E_corr=0;
  PhiIn=0;
  PhiOut=0;

  Dphi=0;
  Theta=0;
  
  m_NumberHigestClusters=0;
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
  m_NumberHigestClusters=m_taggerClusterColl.size();
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
  HigestsClusters.push_back(selectedCluster);
  m_NumberHigestClusters=HigestsClusters.size();
  
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
   if(m_NumberHigestClusters==8)
     {
       // r1 = a1*x + b1
       double a1 = (Y[1]-Y[0])/(X[1]-X[0]);
       double b1 = (Y[0]*X[1]-Y[1]*X[0])/(X[1]-X[0]);
       // r1 = a2*x + b2
       double a2 = (Y[3]-Y[2])/(X[3]-X[2]);
       double b2 = (Y[2]*X[3]-Y[3]*X[2])/(X[3]-X[2]);

       PhiIn   = atan2(Y[1]-Y[0],X[1]-X[0]);
       PhiOut  = atan2(Y[3]-Y[2],X[3]-X[2]);
       Dphi = PhiOut - PhiIn;
       double phiErrIn  = STRIPS_PITCH/(X[1]-X[0]); 
       double phiErrOut = STRIPS_PITCH/(X[3]-X[2]); 
       double DphiErr   = sqrt(phiErrIn*phiErrIn+phiErrOut*phiErrOut);
       
       if(fabs(sin(PhiOut)-sin(PhiIn))>0)
	 {
	   E_rec       = 300.*bt/(sin(PhiOut)-sin(PhiIn));
	   Error_E_rec = 300.*bt/pow(sin(PhiOut)-sin(PhiIn),2.0)*sqrt(pow(cos(PhiOut)*phiErrOut,2.0)+pow(cos(PhiIn)*phiErrIn,2.0));
	 }
       //  Z = A* X + B
       double SX = X[0]+X[1]+X[2]+X[3];
       double SZ = Z[0]+Z[1]+Z[2]+Z[3];
       double SXX = X[0]*X[0] + X[1]*X[1] + X[2]*X[2] + X[3]*X[3];
       double SXZ = X[0]*Z[0] + X[1]*Z[1] + X[2]*Z[2] + X[3]*Z[3];
       double A = (4.*SXZ-SX*SZ)/(4.*SXX-SX*SX);
       double B = (SZ-A*SX)/4.;
       Theta=atan(A);
       if(fabs(cos(Theta))>0)
	 {
	   E_corr = E_rec/cos(Theta);
	   Error_E_corr=Error_E_rec;
	 }
       
       if(fabs(a1-a2)>0)
	 {
	   PX = (b2-b1)/(a1-a2);
	   PY = a1 * PX + b1;
	   PZ = A  * PX + B;
	 }
     }
}

void Recon::report()
{
  std::cout<<" RECON EVENT REPORT: ("<<m_NumberHigestClusters<<")"<<std::endl;
  {
    if(m_NumberHigestClusters==8)
      {
	std::cout<<" Electron Incoming Angle: \t"<<PhiIn<<std::endl;
	std::cout<<" Electron Outgoing Angle: \t"<<PhiOut<<std::endl;
	std::cout<<" Delta angle            : \t"<<Dphi<<std::endl;
	std::cout<<" Theta Angle            : \t"<<Theta<<std::endl;
	std::cout<<" Reconstructed Energy   : \t"<<E_rec<<" +- "<< Error_E_rec <<std::endl;
	std::cout<<" Corrected     Energy   : \t"<<E_corr<<" +- "<< Error_E_corr <<std::endl;
	std::cout<<" Intesection Point      "<<std::endl;
	std::cout<<" \t X \t Y \t Z     "<<std::endl;
	std::cout<<" \t "<<PX<<" \t "<<PY<<" \t "<<PZ<<std::endl;
	std::cout<<" Higest Selected Clusters:    "<<std::endl;
	std::cout<<" \t X \t Y \t Z     "<<std::endl;
	for (int i = 0 ; i < N_MODULES;i++)	   
	  std::cout<<" \t "<<X[i]<<" \t "<<Y[i]<<" \t "<<Z[i]<<std::endl;
	
      }
    else
      {
	std::cout<<" Sorry, not enought clusters!"<<std::endl;
      }
  } 
}

} // end namespace AncillaryData

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
