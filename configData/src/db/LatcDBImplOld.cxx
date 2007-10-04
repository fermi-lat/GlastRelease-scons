
#include "configData/db/LatcDBImplOld.h"
#include <vector>
#include <iostream>


LatcDBImplOld::LatcDBImplOld(){
  m_mq=new MootQuery;
}
const std::string LatcDBImplOld::getFilename(const char* type, unsigned key){
  std::vector<FileDescrip> info;
  std::string types(type);
  std::string rs;
  std::cout<<key<<std::endl;
  std::cout<<"Using old DB implementation"<<std::endl;
  m_mq->getLatcSrc(key,info);
  for (unsigned i=0;i<info.size();i++){
    info[i].getType();
    if (info[i].getType()==types){
      rs= info[i].getPath();
      break;
    }
  }
  //std::cout<<rs<<std::endl;
  return rs;
}
