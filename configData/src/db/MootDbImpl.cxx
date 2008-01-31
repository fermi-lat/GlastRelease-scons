
#include "configData/db/MootDBImpl.h"
#include "mootCore/FileDescrip.h"
#include <vector>
#include <iostream>


MootDBImpl::MootDBImpl(){
  m_mq=new MootQuery;
}

// Local function to generate a list of ParmInfo from an unsigned MOOT key
bool _cfgParmsFromKey(const unsigned configKey, 
                      MootQuery* mq,
                      std::vector<ParmInfo>& cfgParms){
  std::vector<unsigned> pKeys;
  mq->getConfigParmsUsed(configKey, pKeys);
  for ( std::vector<unsigned>::iterator it=pKeys.begin(); it!=pKeys.end(); it++){
    cfgParms.push_back(*mq->getParmInfo(*it));
  }
  return true;
}

// Retrieve a specific parameter file type from a MOOT key
const std::string MootDBImpl::getFilename(const char* type, unsigned key){
  std::vector<ParmInfo> info;
  std::string types(type);
  std::string rs;
  std::cout<<"Using new DB implementation"<<std::endl;
  std::cout<<"Moot key: "<<key<<std::endl;
  _cfgParmsFromKey(key,m_mq,info);
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

// Retrieve the set of used parameter files for a particular MOOT key
const std::vector<std::string> MootDBImpl::getFilenameList(unsigned key){
  std::vector<ParmInfo> info;
  std::string rs;
  std::vector<std::string> rv;
  std::cout<<"Using new DB implementation"<<std::endl;
  std::cout<<"Moot key: "<<key<<std::endl;
  _cfgParmsFromKey(key,m_mq,info);
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
