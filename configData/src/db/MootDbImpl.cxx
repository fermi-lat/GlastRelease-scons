
#include "configData/db/MootDBImpl.h"
#include "configData/fsw/FswEfcSampler.h"
#include "mootCore/FileDescrip.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <set>
#include <cstdlib>
#include "enums/Lsf.h"
#include "facilities/Util.h"


MootDBImpl::MootDBImpl(){
  m_mq=new MootQuery;
  m_acqSummaryCache = 0;
  clearCache();
}


MootDBImpl::~MootDBImpl() {
  delete m_mq;  
  clearCache();
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


MOOT::AcqSummaryInfo* MootDBImpl::getAcqSummary(unsigned startedAt,
						unsigned scid,
						unsigned& configKey ) const {
  if ( startedAt == m_startedAtCache &&
       scid == m_scidCache ) {
    configKey = m_configKeyCache;
    return m_acqSummaryCache;
  }
  clearCache();
  m_acqSummaryCache = m_mq->getAcqSummaryInfo(startedAt,scid);
  if ( m_acqSummaryCache == 0 ) {    
    configKey = LSF_INVALID_UINT;
    return 0;
  }
  try {    
    m_configKeyCache = facilities::Util::stringToUnsigned(m_acqSummaryCache->getConfigKey());
  } catch ( ... ) {
    m_acqSummaryCache = 0;
    m_configKeyCache = LSF_INVALID_UINT;
  }
  return m_acqSummaryCache;
}

unsigned MootDBImpl::prescaleFactor(unsigned configKey,
				    MOOT::LpaMode lpaMode,
				    unsigned handlerId,
				    unsigned filterStatusWord,
				    enums::Lsf::RsdState rsdState,
				    enums::Lsf::LeakedPrescaler leakedPrescaler) const {

  FswEfcSampler* sampler = getFswCdb(configKey,lpaMode,handlerId);
  if ( sampler == 0 ) return LSF_INVALID_UINT;

  switch ( rsdState ) {
  case enums::Lsf::VETOED:
    return 0;
  case enums::Lsf::PASSED:
    return 1;
  case enums::Lsf::IGNORED:
    return sampler->inputPrescaler();
  case enums::Lsf::INVALID:
    return LSF_INVALID_UINT;
  case enums::Lsf::SUPPRESSED:
  case enums::Lsf::LEAKED:
    if ( leakedPrescaler >= 0 ) { return sampler->prescaler(leakedPrescaler); }
    else if ( leakedPrescaler == enums::Lsf::OUTPUT ) { return sampler->outputPrescaler(); }
    else if ( leakedPrescaler == enums::Lsf::INPUT ) { return sampler->inputPrescaler(); }
    break;
  }
  // try and scan the filter status word
  unsigned checkBits = filterStatusWord & sampler->enabled();
  for ( int i(31); i >= 0; i-- ) {
    if (  checkBits & ( 1 << i ) ) return sampler->prescaler(i);
  }
  return 0; 
}


bool MootDBImpl::readEfcFromFile( MOOT::LpaMode lpaMode, unsigned handlerId, const char* fileName ) {
  FswEfcSampler* sampler = FswEfcSampler::makeFromXmlFile(fileName);
  if ( sampler == 0 ) return false;
  if ( lpaMode == LPA_MODE_ALL ) {
    for ( int mode = MOOT::LPA_MODE_NORMAL; mode < MOOT::LPA_MODE_count; mode++ ) {
      unsigned mapKey = FswEfcSampler::mapKey(mode,handlerId);
      m_fswEfcCache[mapKey] = sampler;
    } 
  } else {
    unsigned mapKey = FswEfcSampler::mapKey(lpaMode,handlerId);
    m_fswEfcCache[mapKey] = sampler;
  }
  return true;
}


void MootDBImpl::clearCache() const {
  std::set<FswEfcSampler*> deleted;
  for ( std::map<unsigned,FswEfcSampler*>::iterator itr = m_fswEfcCache.begin(); 
	itr != m_fswEfcCache.end(); itr++ ) {
    if ( deleted.find(itr->second) == deleted.end() ) {
      deleted.insert(itr->second);
      delete itr->second;
    }
  }
  m_fswEfcCache.clear();
  if (m_acqSummaryCache != 0) {
    delete m_acqSummaryCache;
    m_acqSummaryCache = 0;
  }
  m_configKeyCache = LSF_INVALID_UINT;
  m_startedAtCache = LSF_INVALID_UINT;
  m_scidCache = LSF_INVALID_UINT;
}
  
FswEfcSampler* MootDBImpl::getFswCdb(unsigned configKey, MOOT::LpaMode lpaMode, 
				     unsigned handlerId) const {
  if ( configKey != m_configKeyCache ) {
    return 0;
  } 
  unsigned mapKey = FswEfcSampler::mapKey(lpaMode,handlerId);

  std::map<unsigned,FswEfcSampler*>::const_iterator itr = m_fswEfcCache.find(mapKey);
  if ( itr != m_fswEfcCache.end() ) return itr->second;
  std::string handlerName;
  const ConstitInfo* fswConstit = m_mq->getActiveFilter(configKey,lpaMode,handlerId,handlerName);
  if ( fswConstit == 0 ) return 0;

  std::string xmlFilePath = std::string(getenv("MOOT_ARCHIVE"))+"/";
  xmlFilePath += fswConstit->getSrcPath();
  FswEfcSampler* sampler = FswEfcSampler::makeFromXmlFile(xmlFilePath.c_str());
  if ( sampler == 0 ) return 0;
  m_fswEfcCache[mapKey] = sampler;
  return sampler;
}
