#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include "FluxSvc/mainpage.h"
#include "GRBConstants.h"

GRBConstants::GRBConstants(char filen)
{
  ReadParam(filen);
}

void GRBConstants::ReadParam(char filen){
  char buf[100];
  ifstream f1("../src/test/GRBParam.txt");
  if (! f1.is_open()) 
    {
      cout<<" Error Opening Parmas File!!"<<endl;
      cout<<"The file must be placed in:"<<endl;
      cout<<"../src/test/GRBParam.txt"<<endl;
      exit(1);
    }
  cout<<"Read the file: "<<filen<<endl;
  
  f1.getline(buf,100);
  sscanf(buf,"%d",&nshell);
  setNshell(nshell);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&redshift);
  setRedshift(redshift);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&etot);
  setEtot(etot);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&r0);
  setR0(r0);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&t0);
  setT0(t0);
  /*
    fgets(buf,100,f1);
    sscanf(buf,"%lf",&gamma0);
    setGamma0(gamma0);
  */
  f1.close();
}

