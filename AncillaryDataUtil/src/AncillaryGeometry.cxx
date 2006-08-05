#include "AncillaryDataUtil/AncillaryGeometry.h"
#include <iostream>
#include <fstream>

using namespace AncillaryData;

AncillaryGeometry::AncillaryGeometry(std::string GeometryfileName)
{
  FILE *infile;
  infile=fopen(GeometryfileName.c_str(),"r");
  if (infile==NULL)
    {
      std::cerr<<"Error: could not open "<<GeometryfileName<<"\n";
      exit(8);
    }
  fgets(title,sizeof(title),infile);
  //  std::cout<<title<<std::endl;
  fscanf(infile,"%lf ",&BL);
  fgets(title,sizeof(title),infile);
  //std::cout<<title<<std::endl;

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
      std::cout<<v1<<" "<<v2<<" "<<d1<<" "<<d2<<std::endl;
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
}

void AncillaryGeometry::print()
{
  printf("B*L (field times distance in Tm) = %lf\n ",BL);
  std::cout<<title;
  for(unsigned int M=0; M < N_MODULES; M++) 
    std::cout<<M<<"\t "<<getX(M)<<"\t"<<getY(M)<<"\t"<<getZ(M)<<"\t"<<getView1(M)<<"\t"<<getView2(M)<<"\t"<<getDirection1(M)<<"\t"<<getDirection2(M)<<std::endl;
  
}
