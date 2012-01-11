
#include "configData/db/LatcDBImpl.h"
#include "mootCore/FileDescrip.h"
#include <vector>
#include <iostream>
#include <cstdlib>

LatcDBImpl::LatcDBImpl(){
  m_mq=new MootQuery;
}
const std::string LatcDBImpl::getFilename(const char* type, unsigned key){
  std::vector<ParmOffline> info;
  std::string types(type);
  std::string rs;
  std::cout<<"Using new DB implementation"<<std::endl;
  std::cout<<"Moot key: "<<key<<std::endl;
  m_mq->getParmsFromMaster(key,info);
  for (unsigned i=0;i<info.size();i++){
    //std::cout<<info[i].getClass()<<std::endl;
    std::string thepath=info[i].getSrc();
    //std::cout<<thepath<<std::endl;
    if (info[i].getClass()==types){
      rs= thepath;
        break;
    }
  }
  if(rs!="")rs=std::string(getenv("MOOT_ARCHIVE"))+"/"+rs;
  //std::cout<<rs<<std::endl;
  return rs;
}
const std::vector<std::string> LatcDBImpl::getFilenameList( unsigned key){
  std::vector<ParmOffline> info;
  std::string rs;
  std::vector<std::string> rv;
  std::cout<<"Using new DB implementation"<<std::endl;
  std::cout<<"Moot key: "<<key<<std::endl;
  m_mq->getParmsFromMaster(key,info);
  for (unsigned i=0;i<info.size();i++){
    //std::cout<<info[i].getClass()<<std::endl;
    std::string thepath=info[i].getSrc();
    //std::cout<<thepath<<std::endl;
    rs= thepath;
    if(rs!=""){
      rs=std::string(getenv("MOOT_ARCHIVE"))+"/"+rs;
      rv.push_back(rs);
    }
  }
  //std::cout<<rs<<std::endl;
  return rv;
}
