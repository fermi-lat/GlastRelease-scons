// $Header$

// Exercise query routines

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include "mootCore/MoodConnection.h"
#include "mootCore/MootQuery.h"

int main(int nargs, char**)    {

  MOOT::MoodConnection* moodCon = 
    new MOOT::MoodConnection(false);  // read only


  if (!(moodCon) ) {
    std::cerr << std::endl << "unable to connect to MOOD for read "
              << std::endl << std::endl;
    std::cerr.flush();
    return 1;
  }
  else if (!(moodCon->getRdb())) {
    
    std::cerr << std::endl << "unable to connect to MOOD for read  "
              << std::endl << std::endl;
    std::cerr.flush();
    return 1;
  }

  else {
    std::cout << std::endl << "Successfully connected "
              << " for read " << std::endl << std::endl 
              << std::endl << std::endl;
  }
  
  MOOT::MootQuery q(moodCon);
  
  std::vector<std::string> inputKeys;
  //  bool gottem = q.getConfigInputsUn(4, inputKeys);
  bool gottem = q.getConfigInputs(4, inputKeys);
    
  if (!gottem) {
    std::cout << "MootQuery::getConfigInputs  failed" << std::endl;
  }
  else {
    std::cout << "MootQuery::getConfigInputs for config #4 returns: "
              << std::endl;
    for (unsigned iVal = 0; iVal < inputKeys.size(); iVal++) {
      std::cout << inputKeys[iVal] << std::endl;
    }
    std::cout << std::endl;
    std::cout.flush();
  }

  std::vector<std::string> paths;
  gottem = q.getConfigFmxPaths(4, paths);

  
  if (!gottem) {
    std::cout << "MootQuery::getConfigFmxPaths failed" << std::endl;
  }
  else {
    std::cout << "MootQuery::getConfigFmxPaths for config #4 returns: "
              << std::endl;
    for (unsigned iVal = 0; iVal < paths.size(); iVal++) {
      std::cout << paths[iVal] << std::endl;
    }
    std::cout << std::endl;
    std::cout.flush();
  }

  paths.clear();
  gottem = q.getConfigFmxPathsByName("configOf2006-05-26 23:04:16", paths);
  
  if (!gottem) {
    std::cout << "MootQuery::getConfigFmxPathsByName  failed" << std::endl;
  }
  else {
    std::cout << "MootQuery::getConfigFmxPathsByName for config " 
              << "'configOf2006-05-26 23:04:16' returns: "
              << std::endl;
    for (unsigned iVal = 0; iVal < paths.size(); iVal++) {
      std::cout << paths[iVal] << std::endl;
    }
    std::cout << std::endl;
    std::cout.flush();
  }

  paths.clear();
  gottem = q.getConfigFmxPathsByAlg("none", "1", paths);
  
  if (!gottem) {
    std::cout << "MootQuery::getConfigFmxPathsByAlg  failed" << std::endl;
  }
  else {
    std::cout << "MootQuery::getConfigFmxPathsByAlg for algorithm " 
              << "'none' returns: "
              << std::endl;
    for (unsigned iVal = 0; iVal < paths.size(); iVal++) {
      std::cout << paths[iVal] << std::endl;
    }
    std::cout << std::endl;
    std::cout.flush();
  }

  std::vector<MOOT::FileDescrip> descrips;
  gottem = q.getLatcSrc("229", descrips);
  if (!gottem) {
    std::cout << "MootQuery::getLatcSrc failed " << std::endl;
  }
  else {
    std::cout << "MootQuery::getLatcSrc for fmx key 229 returns: " 
              << std::endl;
    for (unsigned iFile = 0; iFile < descrips.size(); iFile++) {
      std::cout << "File " << iFile << ":" << std::endl <<
        descrips[iFile].m_path << std::endl 
                << " of type " << descrips[iFile].m_type << std::endl;
    }
    std::cout << std::endl;
    std::cout.flush();
  }

  std::vector<MOOT::FswInfo> info;
  gottem = q.getConfigFswInfo("4", info);
  if (!gottem) {
    std::cout << "MootQuery::getConfigFswInfo failed " << std::endl;
  }
  else {
    std::cout << "MootQuery::getFswInfo for config key 4 returns: " 
              << std::endl;
    for (unsigned i = 0; i < info.size(); i++) {
      std::cout << "fsw input " << i << ":" << std::endl
                << "fmx path: " << info[i].getFmxPath() << std::endl 
                << "src: " << info[i].getSrc() << std::endl  
                << "description: " << info[i].getDescription() << std::endl
                << "fmx key: " << info[i].getFmxKey() << std::endl 
                << "MOOT class key: " << info[i].getMootClassKey() 
                << std::endl  << std::endl;
    }
    std::cout << std::endl;
    std::cout.flush();
  }

  std::cout << "***** New queries (8 Sept. 2006) *******" << std::endl;
  std::vector<unsigned> keys;

  std::string cname("secondConfig");

  std::string algname("emptyAlg");
  std::string dontcare("");

  unsigned nKeys = q.listConfigKeys(keys);

  std::cout << "Config keys for configs with status = 'CREATED': ";
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;

  keys.clear();
  nKeys = q.listConfigKeys(keys, "STARTED");

  std::cout << "Config keys for configs with status = 'STARTED': ";
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;

  keys.clear();
  nKeys = q.listConfigKeys(keys, dontcare);

  std::cout << "All config keys: ";
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;
  
  unsigned theKey = q.getLastConfigKeyByName(cname);
  std::cout << "Last good key with name '" << cname << "': "
            << theKey << std::endl;

  theKey = q.getLastConfigKeyByName(cname, dontcare);
  std::cout << "Last  key with name '" << cname << "': "
            << theKey << std::endl;

  theKey = q.getLastConfigKeyByAlg(algname);
  std::cout << "Last good key with alg name '" << algname << "', step 1: "
            << theKey << std::endl;

  theKey = q.getLastConfigKeyByAlg(algname, 1, std::string("INVALID"));
  std::cout << "Last INVALID key with alg name '" << algname << "', step 1: "
            << theKey << std::endl;


  theKey = q.getLastConfigKeyByAlg(algname, 0, dontcare);
  std::cout << "Last key (good, bad or indifferent) with alg name '" 
            << algname << "', any step: "
            << theKey << std::endl;

  keys.clear();
  nKeys = q.getConfigKeysByName(keys, cname);
  std::cout << "Config keys for configs with status = 'CREATED', name "
            << cname;
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;

  keys.clear();
  nKeys = q.getConfigKeysByName(keys, cname, dontcare);
  std::cout << "Config keys for configs with any status, name "
            << cname;
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;

  keys.clear();
  nKeys = q.getConfigKeysByAlg(keys, algname);
  std::cout << "Config keys for configs with CREATED status, alg name "
            << algname << ", first step";
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;


  keys.clear();
  nKeys = q.getConfigKeysByAlg(keys, algname, 2, dontcare);
  std::cout << "Config keys for configs with any status, alg name "
            << algname << ", 2nd step";
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;

  keys.clear();
  nKeys = q.getConfigKeysByAlg(keys, algname, 0, dontcare, dontcare);
  std::cout << "Config keys for configs with any status, any instrument,"
            << " alg name " << algname << ", any step";
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;


  keys.clear();
  nKeys = q.getConfigKeysByAlg(keys, algname, 0, dontcare, 
                               std::string("violin"));
  std::cout << "Config keys for configs with any status, instr='violin',"
            << " alg name " << algname << ", any step";
  for (unsigned i = 0; i < nKeys; i++) {
    if ((i % 10) == 0) std::cout << std::endl;
    std::cout << keys[i] << " ";
  }
  std::cout << std::endl;

  std::vector<MOOT::ConfigInfo> cInfo;

  unsigned nConfigs = q.getConfigInfo(cInfo);
  MOOT::ConfigInfo* pInfo = &cInfo[0];

  std::cout << "Retrieved information for " << nConfigs << " Configs "
            << std::endl;
  std::cout << "For first config: " << std::endl;
  std::cout << "Key = " << pInfo->getKey() << std::endl;
  std::cout << "Name = " << pInfo->getName() << std::endl;
  std::cout << "Alg. = " << pInfo->getAlg() << std::endl;
  std::cout << "Step = " << pInfo->getStep() << std::endl;
  std::cout << "Description = " << pInfo->getDescrip() << std::endl;
  std::cout << "Status = " << pInfo->getStatus() << std::endl;
  std::cout << "active_state = " << pInfo->getActive() << std::endl;
  std::cout << "Mode = " << pInfo->getMode() << std::endl;
  std::cout << "Creation time = " << pInfo->getCreationTime() 
            << std::endl << std::endl;

  pInfo = &cInfo[nConfigs - 1];

  std::cout << "For last config: " << std::endl;
  std::cout << "Key = " << pInfo->getKey() << std::endl;
  std::cout << "Name = " << pInfo->getName() << std::endl;
  std::cout << "Alg. = " << pInfo->getAlg() << std::endl;
  std::cout << "Step = " << pInfo->getStep() << std::endl;
  std::cout << "Description = " << pInfo->getDescrip() << std::endl;
  std::cout << "Status = " << pInfo->getStatus() << std::endl;
  std::cout << "Mode = " << pInfo->getMode() << std::endl;
  std::cout << "Creation time = " << pInfo->getCreationTime() 
            << std::endl  << std::endl;

  cInfo.clear();
  nConfigs = q.getConfigInfo(cInfo, "", "");
  std::cout << "Using no cuts at all, found " << nConfigs << " configs" 
            << std::endl;
  nConfigs = q.getConfigInfo(cInfo, "STARTED", "");  
 
  std::cout << "Found " << nConfigs << " with status = STARTED" << std::endl;
  return 0;
}
