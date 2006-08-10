#include "AncillaryDataEvent/Recon.h"
#include <algorithm>

//using namespace AncillaryData;

namespace AncillaryData {

Recon::Recon(AncillaryData::Digi *digiEvent)
{
  setEventNumber(digiEvent->getEventNumber());
  setSpillNumber(digiEvent->getSpillNumber());
  setQdcHitCol(digiEvent->getQdcHitCol());      
  setScalerHitCol(digiEvent->getScalerHitCol());      

  PX = -9999.0;  PY= -9999.0;  PZ= -9999.0;
  E_rec=0; 
  E_corr=0;
  PhiIn=-100;
  PhiOut=-100;

  Dphi=-100;
  Theta=-100;
  
  m_NumberHigestClusters=0;
  m_NumberTotalClusters=0;
  for(unsigned int m=0; m < N_MODULES; m++)
    {
      m_NumberClusters[0][m]=0;
      m_NumberClusters[1][m]=0;
    }


}

void Recon::print()
{
  std::cout<< " Ancillary Recon Event: "<<getEventNumber()<<" Spill Number: "<<getSpillNumber()<<std::endl;
  std::cout<< " --- number of Tagger Clusters: "<<m_taggerClusterCol.size()<<std::endl;
  for(std::vector<TaggerCluster>::iterator pos=m_taggerClusterCol.begin(); pos!=m_taggerClusterCol.end(); ++pos)
    (*pos).print();
  std::cout<< " --- number of QDC Hits       : "<<m_qdcHitCol.size()<<std::endl;
  for(std::vector<QdcHit>::iterator        pos=       m_qdcHitCol.begin(); pos!=       m_qdcHitCol.end(); ++pos)
    (*pos).print();
  std::cout<< " --- number of Scaler Hits       : "<<m_scalerHitCol.size()<<std::endl;
  for(std::vector<ScalerHit>::iterator        pos=       m_scalerHitCol.begin(); pos!=       m_scalerHitCol.end(); ++pos)
    (*pos).print();

}


std::vector<TaggerCluster> Recon::GetHighestClusters()
{
  SortClusters();
  for(unsigned int m=0;m < N_MODULES;m++)
    {
      m_NumberClusters[0][m]=0;
      m_NumberClusters[1][m]=0;
    }
  //  std::cout<<"std::vector<TaggerCluster> Recon::GetHighestClusters: "<<m_taggerClusterCol.size()<<std::endl;
  //  m_NumberClusters       = m_taggerClusterCol.size();
  if(m_taggerClusterCol.size()<=1) 
    return m_taggerClusterCol;
  std::vector<TaggerCluster> HigestsClusters;
  ///////
  for(std::vector<TaggerCluster>::iterator pos=m_taggerClusterCol.begin(); pos<m_taggerClusterCol.end();pos++)
    m_NumberClusters[(*pos).getLayerId()][(*pos).getModuleId()]++;
  m_NumberTotalClusters = 0;
  for(unsigned int m=0; m < N_MODULES; m++)
    {
      m_NumberTotalClusters += m_NumberClusters[0][m];
      m_NumberTotalClusters += m_NumberClusters[1][m];
    }
  //////
  std::vector<TaggerCluster>::iterator pos=m_taggerClusterCol.begin();
  TaggerCluster selectedCluster=(*pos);
  pos++;
  while(pos<m_taggerClusterCol.end())
    {
      TaggerCluster newCluster=(*pos);
      //      std::cout<<"selected cluster:"<<selectedCluster.getModuleId()<<", "<<selectedCluster.getLayerId()<<std::endl;
      //      std::cout<<"new Cluster cluster:"<<newCluster.getModuleId()<<", "<<newCluster.getLayerId()<<std::endl;

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
  //  std::cout<<"std::vector<TaggerCluster> Recon::GetHighestClusters m_NumberHigestClusters "<<m_NumberHigestClusters<<std::endl;
  return HigestsClusters;
}

void Recon::SortClusters()
{
  std::sort(m_taggerClusterCol.begin(),m_taggerClusterCol.end(),ClusterSortPredicate);
}


void Recon::computePositions(AncillaryGeometry *geometry)
{
  std::vector<TaggerCluster> higestClusters = GetHighestClusters();
  for (unsigned int i = 0 ; i < N_MODULES ; i++)
    {
      X[i] = 0.0; Y[i] = 0.0; Z[i] = 0.0;
    }
  
  for(std::vector<TaggerCluster>::iterator pos=higestClusters.begin();pos!=higestClusters.end(); ++pos)
    {
      const unsigned int M= (*pos).getModuleId();
      const unsigned int L= (*pos).getLayerId(); // 0 or 1
      // View and direction of the first view
      unsigned int V,D;
      if(L==0)
	{
	  V = geometry->getView1(M);
	  D = geometry->getDirection1(M);
	}      
      else
	{
	  V = geometry->getView2(M);
	  D = geometry->getDirection2(M);
	}
      const double HalfWafer = N_CHANNELS_PER_LAYER*STRIPS_PITCH/2.0;
      // this is valid for both views:
      X[M] = geometry->getX(M);
      // case of first view
      if(V==1 && D == 0)      // case : V = Z, D = +
	Z[M] = geometry->getZ(M) - HalfWafer + (*pos).getPosition();
      else if(V==1 && D == 1 )      // case : V = Z, D = -
	Z[M] = geometry->getZ(M) + HalfWafer - (*pos).getPosition();
      else if(V==0 && D == 0)      // case : V = Y, D = +
	Y[M] = geometry->getY(M) - HalfWafer + (*pos).getPosition();
      else if(V==0 && D == 1 )      // case : V = Y, D = -
	Y[M] = geometry->getY(M) + HalfWafer - (*pos).getPosition();
    }
}

void Recon::reconstructEnergy(AncillaryGeometry *geometry)
{
   E_rec  = 0.0;
   E_corr = 0.0;
   double bt=geometry->getBL();
   double beamMomentum=geometry->getBeamMomentum()*1000.0; //MeV
   
   // first track:
   const double Dist1 = X[1]-X[0];
   const double Disp1 = Y[1]-Y[0];
   const double Dist2 = X[3]-X[2];
   const double Disp2 = Y[3]-Y[2];
   double a1,a2,b1,b2;

   if(Disp1!=0)
     {
       a1 = Disp1/Dist1;
       b1 = (Y[0]*X[1]-Y[1]*X[0])/Dist1;
       PhiIn   = atan2(Disp1,Dist1);
     }
   
   if(Disp2!=0)
     {
       // r1 = a2*x + b2
       a2 = Disp2/Dist2;
       b2 = (Y[2]*X[3]-Y[3]*X[2])/Dist2;
       PhiOut  = atan2(Disp2,Dist2);
     }
   
   if(Disp1!=0 && Disp2!=0 )
     {
       Dphi    = PhiOut - PhiIn;
       double phiErrIn  = STRIPS_PITCH/Dist1;
       double phiErrOut = STRIPS_PITCH/Dist2;
       //       double DphiErr   = sqrt(phiErrIn*phiErrIn+phiErrOut*phiErrOut);
       if(fabs(sin(PhiOut)-sin(PhiIn))>0)
	 {
	   E_rec       = beamMomentum - (300.*bt/(sin(PhiOut)-sin(PhiIn)));
	   Error_E_rec = (300.*bt/pow(sin(PhiOut)-sin(PhiIn),2.0)*sqrt(pow(cos(PhiOut)*phiErrOut,2.0)+pow(cos(PhiIn)*phiErrIn,2.0)));
	   if(fabs(a1-a2)>0)
             {
               PX = (b2-b1)/(a1-a2);
               PY = a1 * PX + b1;
	     }   
	 }
       //  Z = A* X + B
       // Case with two points in z:
       
       if(m_NumberHigestClusters==8 && m_NumberTotalClusters == 8)
	 {
	   double SX = X[0]+X[1]+X[2]+X[3];
	   double SZ = Z[0]+Z[1]+Z[2]+Z[3];
	   double SXX = X[0]*X[0] + X[1]*X[1] + X[2]*X[2] + X[3]*X[3];
	   double SZX = Z[0]*X[0] + Z[1]*X[1] + Z[2]*X[2] + Z[3]*X[3];
	   double A = (4.*SZX-SZ*SX)/(4.*SXX-SX*SX);
	   double B = (SZ-A*SX)/4.;
	   if(fabs(a1-a2)>0) PZ = A  * PX + B;
	   Theta=atan2(4.*SZX-SZ*SX,4.*SXX-SX*SX);
	   if(fabs(cos(Theta))>0)
	     {
	       E_corr = E_rec/cos(Theta);
	       Error_E_corr=Error_E_rec;
	     }
	 }
     }
}

void Recon::report()
{
  std::cout<<" RECON EVENT REPORT: Total Clusters: "<<m_NumberTotalClusters<<" Higest:" <<m_NumberHigestClusters<<std::endl;
  {
    const double Disp1 = Y[1]-Y[0];
    const double Disp2 = Y[3]-Y[2];
    
    if(Disp1!=0)
      std::cout<<" Electron Incoming Angle: \t"<<PhiIn<<std::endl;
    if(Disp2!=0)
      std::cout<<" Electron Outgoing Angle: \t"<<PhiOut<<std::endl;
    if(Disp1!=0 && Disp2!=0)
      {
	std::cout<<" Delta angle            : \t"<<Dphi<<std::endl;
	std::cout<<" Reconstructed Energy   : \t"<<E_rec<<" +- "<< Error_E_rec <<std::endl;
	if(m_NumberHigestClusters==8 && m_NumberTotalClusters == 8)
	  {
	    std::cout<<" Theta Angle            : \t"<<Theta<<std::endl;
	    std::cout<<" Corrected     Energy   : \t"<<E_corr<<" +- "<< Error_E_corr <<std::endl;
	  }
	std::cout<<" Intesection Point      "<<std::endl;
	std::cout<<" \t X \t Y \t Z     "<<std::endl;
	std::cout<<" \t "<<PX<<" \t "<<PY<<" \t "<<PZ<<std::endl;
	std::cout<<" Higest Selected Clusters:    "<<std::endl;
      }
    std::cout<<" \t X \t Y \t Z     "<<std::endl;
    for(unsigned int m=0; m < N_MODULES; m++)
      std::cout<<" \t "<<X[m]<<" \t "<<Y[m]<<" \t "<<Z[m]<<std::endl;
    /*
      std::cout<<" Sorry, not enough clusters!"<<std::endl;
      std::cout<<" M = "<<m<<" Num Clusters: L = 0: "<< m_NumberClusters[0][m]<<" L = 1: "<<m_NumberClusters[1][m]<<std::endl;
    */
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
  TAG_EGAMMA_ERR 	 error on photon energy measuement
  CRNKV_PHA[2] 	 pulse height amplitude in the 2 cerenkov
  SCINT_PHA[8] 	 array of floats - PHA for scintillators in the setup (8 is a placeholder for a reasonable number, might change)
  TRD_NUM[16]

  void Recon::computeFinalTupla()
  {
  
  }
*/
