#include "AncillaryDataUtil/AncillaryGeometry.h"
#include <iostream>
#include <fstream>

using namespace AncillaryData;

AncillaryGeometry::AncillaryGeometry(std::string GeometryfileName)
{
  char FILENAME[100];
  sprintf(FILENAME,GeometryfileName.c_str());
  FILE *infile;
  infile=fopen(FILENAME,"r");
  if (infile==NULL)
    {
      std::cerr<<"Error: could not open "<<FILENAME<<"\n";
      exit(8);
    }
  fgets(title,sizeof(title),infile);
  //  std::cout<<title<<std::endl;
  fscanf(infile,"%lf ",&BL);
  fgets(title,sizeof(title),infile);
  //std::cout<<title<<std::endl;

  unsigned int id;
  double x,y,z,tx,ty,tz;
  char v[1];
  unsigned int M=0;
  for(unsigned int i=0; i < N_MODULES; i++) 
    {
      fscanf(infile,"%u %lf %lf %lf %lf %lf %lf %c",&id,&x,&y,&z,&tx,&ty,&tz,&v); 
      M=id;
      X[M]=x;
      Y[M]=y;
      Z[M]=z;
      Tx[M]=tx;
      Ty[M]=ty;
      Tz[M]=tz;
      if(v=="Y" || "y")  View[M]=0; // not flipped (strip along Z)
      else if(v=="Z" || "z")  View[M]=0; // flipped (strip along Y)
      else std::cout<<"Wrong geometry file format"<<std::endl;
    }
}

void AncillaryGeometry::print()
{
  printf("B*L (field times distance in Tm) = %lf",BL);
  std::cout<<title;
  for(unsigned int M=0; M < N_MODULES; M++) 
    std::cout<<M<<"\t "<<X[M]<<"\t"<<Y[M]<<"\t"<<Z[M]<<"\t"<<Tx[M]<<"\t"<<Ty[M]<<"\t"<<Tz[M]<<"\t"<<View[M]<<std::endl;
  
}
