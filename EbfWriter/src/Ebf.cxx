#include "EbfWriter/Ebf.h"

using namespace EbfWriterTds;

Ebf::Ebf(){
  data=NULL;
  length=0;
}

Ebf::Ebf(char *newData, unsigned int dataLength){
  data=NULL;
  length=0;
  set(newData,dataLength);
}

Ebf::~Ebf(){
  if(data!=NULL)
    delete[] data;
}

char *Ebf::get(unsigned int &dataLength) const{
  dataLength=length;
  return data;
}

inline std::ostream& Ebf::fillStream( std::ostream &s) const{
  if(length>0){
    s.write(data,length);
  }
  return s;
}

std::ostream& operator << (std::ostream& s, const Ebf& obj){
  return obj.fillStream(s);
}

void Ebf::set(char *newData,unsigned int dataLength){
  if(data!=NULL)
    delete[] data;
  data=NULL;
  data=new char[dataLength];
  memcpy(data,newData,dataLength);
  length=dataLength;
}
