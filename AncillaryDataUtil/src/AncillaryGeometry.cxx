#include "AncillaryDataUtil/AncillaryGeometry.h"
#include <iostream>
#include <fstream>
#include "facilities/Util.h"

using namespace AncillaryData;

AncillaryGeometry::AncillaryGeometry(std::string GeometryfileName,std::string rcReportPath)
{
  facilities::Util::expandEnvVar(&GeometryfileName);
  facilities::Util::expandEnvVar(&rcReportPath);
  FILE *infile;
  infile=fopen(GeometryfileName.c_str(),"r");
  if (infile==NULL)
    {
      std::cerr<<"Error: could not open "<<GeometryfileName<<"\n";
      exit(8);
    }
  fgets(title,sizeof(title),infile);
  fscanf(infile,"%lf %lf",&BL0,&BL1);
  char tmp[200];
  fgets(tmp,sizeof(tmp),infile);
  fgets(tmp,sizeof(tmp),infile);

  unsigned int id;
  double x,y,z;
  char v1,v2,d1,d2;
  unsigned int M=0;
  for(unsigned int i=0; i < N_MODULES; i++) 
    {
      fscanf(infile,"%u %lf %lf %lf %c %c %c %c",&id,&x,&y,&z,&v1,&d1,&v2,&d2);
      M=id;
      X[M]=x;
      Y[M]=y;
      Z[M]=z;
      //      std::cout<<id<<" "<<x<<" "<<y<<" "<<z<<" "<<v1<<" "<<v2<<" "<<d1<<" "<<d2<<std::endl;
      if (v1=='Y'|| v1=='y') V1[M]=0;
      else if (v1=='Z'|| v1=='z') V1[M]=1;
      else std::cout<<"Wrong geometry file format"<<std::endl;
      if (v2=='Y'|| v2=='y') V2[M]=0;
      else if (v2=='Z'|| v2=='z') V2[M]=1;
      else std::cout<<"Wrong geometry file format"<<std::endl;
      if (d1=='+') D1[M]=0;
      else if (d1=='-') D1[M]=1;
      else std::cout<<"Wrong geometry file format"<<std::endl;
      if (d2=='+') D2[M]=0;
      else if (d2=='-') D2[M]=1;
      else std::cout<<"Wrong geometry file format"<<std::endl;
    }
  std::ifstream rcReport(rcReportPath.c_str());
  if(!rcReport.is_open())
    {
      m_MagnetCurrent=1.0;
      BL=1.0;
      m_BeamMomentum=0.0;
      std::cerr<<"WARNING: AncillaryUtil:: could not open "<<rcReportPath<<"\n";
    }
  else
    {
      std::string input;
      getline(rcReport, input, '\n');
      int pos1=input.rfind("<MagnetCurrent>");
      int pos2=input.rfind("</MagnetCurrent>");
      std::string magnetic_current_value = input.substr(pos1+15, pos2-pos1);
      pos1=input.rfind("<BeamMomentum>");
      pos2=input.rfind("</BeamMomentum>");
      std::string beam_momuentum = input.substr(pos1+14, pos2-pos1);
      m_MagnetCurrent = atof(magnetic_current_value.c_str());
      m_BeamMomentum = atof(beam_momuentum.c_str());
    }
  BL = ((m_MagnetCurrent==0) ? BL = 1.0 : BL0 * m_MagnetCurrent + BL1 * m_MagnetCurrent*m_MagnetCurrent);
}

void AncillaryGeometry::print()
{
  std::cout<<title<<std::endl;
  for(unsigned int M=0; M < N_MODULES; M++) 
    std::cout<<M<<"\t "<<getX(M)<<"\t"<<getY(M)<<"\t"<<getZ(M)<<"\t"<<getView1(M)<<"\t"<<getView2(M)<<"\t"<<getDirection1(M)<<"\t"<<getDirection2(M)<<std::endl;
  std::cout<<"------------------------------------------------------------------------"<<std::endl;
  if(m_MagnetCurrent==0) 
    {
      std::cout<<" Magent Current =0; BL scaled to 1, Beam Momentum: "<<m_BeamMomentum<<" GeV"<<std::endl;
    }
  else 
    {
      std::cout<<"B*L (field times distance in Tm) =  "<<BL0<<" * "<<m_MagnetCurrent<<" + "<<BL1<<" *("<<m_MagnetCurrent<<")^2 = "<<BL<<std::endl;
      std::cout<<" Beam Momentum: "<<m_BeamMomentum<<" GeV"<<std::endl;
    }
}

