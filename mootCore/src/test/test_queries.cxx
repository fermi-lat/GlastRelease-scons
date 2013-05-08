// $Header$

// Exercise query routines

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include "mootCore/MoodConnection.h"
#include "mootCore/MootQuery.h"
#include "facilities/Timestamp.h"
#include "facilities/commonUtilities.h"

void writeInfo(MOOT::ConfigInfo* pInfo, std::ostream& out);
void writeConstituent(const MOOT::ConstitInfo* pInfo, std::ostream& out);
void writeParmOff(const std::vector<MOOT::ParmOffline>& parm);
void writeAcqInfo(const MOOT::AcqSummaryInfo* pInfo, std::ostream& out);

int main(int /* nargs */, char**)    {

  facilities::commonUtilities::setupEnvironment();

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
  //  bool gottem;
  /*
  try {
    gottem = q.getConfigInputs(4, inputKeys);
  }
  catch (std::exception ex) {
    std::cerr << "failed with exception" << std::endl;
    gottem = false;
  }
  catch (MOOT::DbUtilException mex) {
    std::cerr << "failed with MOOT::DbUtilException" << std::endl;
    gottem = false;
  }
    
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
  */
  std::vector<std::string> paths;
  /*
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
  */
  /*
  std::cout << "***** New queries (8 Sept. 2006) *******" << std::endl;
  */

  std::vector<unsigned> keys;
  // was commented out
  /* and is again
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
  */
  // end of "was commented out" section
  std::vector<MOOT::ConfigInfo> cInfo;

  unsigned nConfigs = q.getConfigInfo(cInfo);
  MOOT::ConfigInfo* pInfo = &cInfo[0];

  std::cout << "Retrieved information for " << nConfigs << " Configs "
            << std::endl;
  std::cout << "For first config: " << std::endl;

  writeInfo(pInfo, std::cout);

  pInfo = &cInfo[nConfigs - 1];

  std::cout << "For last config: " << std::endl;
  writeInfo(pInfo, std::cout);

  cInfo.clear();
  nConfigs = q.getConfigInfo(cInfo, "", "");
  std::cout << "Using no cuts at all, found " << nConfigs << " configs" 
            << std::endl;
  nConfigs = q.getConfigInfo(cInfo, "STARTED", "");  
 
  std::cout << "Found " << nConfigs << " with status = STARTED" << std::endl;
  
  // Exercise getConfigsRequested, getConfigsUsed
  /*
  std::vector<unsigned> parmKeys;
  parmKeys.reserve(10);

  unsigned configKey = 119;
  bool ok;

  for (unsigned itry = 0; itry < 2; itry++) {
    try {
      ok = q.getConfigParmsRequest(configKey, parmKeys);
    }
    catch (std::exception ex) {
      std::cerr << "getConfigParmsRequest of key #" << configKey
                << " failed " << std::endl;
      ok = false;
    }
    if (ok) {
      std::cout << "Config #" << configKey << " requested parameters ";
      for (unsigned i = 0; i < parmKeys.size(); i++) {
        if ((i % 10) == 0) std::cout << std::endl;
        std::cout << parmKeys[i] << " ";
      }
      std::cout << std::endl << std::endl;
    }
    
    parmKeys.clear();

    try {
      ok = q.getConfigParmsUsed(configKey, parmKeys);
    }
    catch (std::exception ex) {
      std::cerr << "getConfigParmsUsed of key #" << configKey
                << " failed with error " << ex.what() << std::endl;
      ok = false;
    }
    if (ok) {
      std::cout << "Config #" << configKey << " used parameters ";
      for (unsigned i = 0; i < parmKeys.size(); i++) {
        if ((i % 10) == 0) std::cout << std::endl;
        std::cout << parmKeys[i] << " ";
      }
      std::cout << std::endl << std::endl;
    }
    configKey = 120;
    parmKeys.clear();
  }

  unsigned ckey = 8;
  pInfo = q.getConfigInfo(ckey);
  if (pInfo) {
    std::cout << "Got config info for config #" << ckey << std::endl;
    writeInfo(pInfo, std::cout);
  }
  else {
    std::cout << "No info found for config #" << ckey 
              << std::endl << std::endl;
  }

  ckey = 119;
  pInfo = q.getConfigInfo(ckey);
  if (pInfo) {
    std::cout << "Got config info for config #" << ckey << std::endl;
    writeInfo(pInfo, std::cout);
  }
  else {
    std::cout << "No info found for config #" << ckey 
              << std::endl << std::endl;
  }

  ckey = 10023;
  pInfo = q.getConfigInfo(ckey);
  if (pInfo) {
    std::cout << "Got config info for config #" << ckey << std::endl;
    writeInfo(pInfo, std::cout);
  }
  else {
    std::cout << "No info found for config #" << ckey 
              << std::endl << std::endl;
  }
  */
  std::cout << "***** New queries (10 May 2007) *******" << std::endl;

  keys.clear();
  q.listAncAliasKeys(keys);

  std::cout << "listAncAliasKeys(keys) returned " << keys.size() 
            << " key(s)" << std::endl;

  keys.clear();
  q.listAncAliasKeys(keys, "ACD_PedestalCalibration");

  std::cout << 
    "listAncAliasKeys(keys,'ACD_PedestalCalibration') returned " 
            << keys.size() 
            << " key(s)" << std::endl;

  keys.clear();
  q.listAncAliasKeys(keys, "ACD_HldCalibration");

  std::cout << 
    "listAncAliasKeys(keys,'ACD_HldCalibration') returned " 
            << keys.size() 
            << " key(s)" << std::endl;

  keys.clear();
  try {
    q.listAncAliasKeys(keys, "bad_class_name");
    std::cout << 
      "listAncAliasKeys(keys,'bad_class_name') returned " 
              << keys.size() 
              << " key(s)" << std::endl;
  }
  catch (std::exception ex) {
    std::cerr << 
      "listAncAliasKeys(keys,'bad_class_name') threw an exception " 
              << std::endl;
  }
  /*
  keys.clear();
  q.listVoteKeys(keys);

  std::cout << "listVote(keys) returned " << keys.size() 
            << " key(s)" << std::endl;

  keys.clear();
  q.listVoteKeys(keys, "generic");

  std::cout << "listVote(keys, 'generic') returned " << keys.size() 
            << "keys" << std::endl;

  keys.clear();
  q.listVoteKeys(keys, "ACD_PHA");

  std::cout << "listVote(keys, 'ACD_PHA') returned " << keys.size() 
            << "keys" << std::endl;

  keys.clear();
  q.listVoteKeys(keys, "ACD_PHA", "*");

  std::cout << "listVote(keys, 'ACD_PHA', '*') returned " << keys.size() 
            << " key(s)" << std::endl;

  keys.clear();
  q.listVoteKeys(keys, "ACD_PHA", "INVALID");

  std::cout << "listVote(keys, 'ACD_PHA', 'INVALID') returned " << keys.size() 
            << " key(s)" << std::endl;
  */
  std::vector<std::string> names;

  names.clear();
  q.getPrecincts(names); 

  std::cout << names.size() << " precinct names found.  They are:" 
            << std::endl;
  for (unsigned ix = 0; ix < names.size(); ix++) {
    std::cout << names[ix] << " ";
    if ((ix % 6) == 5) std::cout << std::endl;
  }
  std::cout << std::endl;

  names.clear();
  q.getAncClasses(names); 

  std::cout << names.size() << " anc. class names found.  They are:" << std::endl;
  for (unsigned ix = 0; ix < names.size(); ix++) {
    std::cout << names[ix] << " ";
    if ((ix % 3) == 2) std::cout << std::endl;
  }
  std::cout << std::endl;

  names.clear();
  std::string prec("ACD_Mode");
  q.getParmClasses(names, prec);

  std::cout << names.size() << " parm classes names belonging to precinct '"
            << prec <<  "' found.  They are:" << std::endl;
  for (unsigned ix = 0; ix < names.size(); ix++) {
    std::cout << names[ix] << " ";
    if ((ix % 3) == 2) std::cout << std::endl;
  }
  std::cout << std::endl;


  std::vector<unsigned> pkeys;
  unsigned voteKey;
  //  unsigned minKey = 1;
  unsigned minKey = 16;
  unsigned maxKey = 18;
  for (voteKey = minKey; voteKey < maxKey; voteKey++) {
    bool isUpToDate = q.getVoteParameters(voteKey, pkeys);
    if (isUpToDate) {
      std::cout << "For vote with key=" << voteKey << " found " << pkeys.size()
                << " parameter files with keys:" << std::endl;
      for (unsigned ikey = 0; ikey < pkeys.size(); ikey++) {
        std::cout << pkeys[ikey] << std::endl;
      }
    }
    else std::cout << "Vote with key=" << voteKey << " was not up to date" 
                   << std::endl;
  }

  std::vector<MOOT::ParmOffline> parmOff;
  unsigned fmxMasterId = 4617;

  bool done = true;
  do {
    bool getStatus = q.getParmsFromMaster(fmxMasterId, parmOff);
    done = true;

    if (!getStatus) {
      std::cout << "MootQuery::getParmsFromMaster failed for fmx logical id = "    
                << fmxMasterId << std::endl;
    }
    else  {
      std::cout << "MootQuery::getParmsFromMaster for fmx logical id = "    
                << fmxMasterId << " returns info for " << parmOff.size()
                << " parameter files:" << std::endl;
      writeParmOff(parmOff);
    }
  } while (!done) ;

  /*   Not supported on prod. dbs yet
  MOOT::ConstitInfo* byKey = 
    q.getConstituentInfo(75);
  writeConstituent(byKey, std::cout);

  MOOT::ConstitInfo* byId = q.getConstituentByFswId(6025);
  writeConstituent(byId, std::cout);

  
  std::cout << std::endl << "Output from getLpaConstituents: " << std::endl;
  std::vector<MOOT::ConstitInfo> lpa;
  bool ok = q.getLpaConstituents(6241, lpa);
  if (ok) {
    for (unsigned ix = 0; ix < lpa.size(); ix++) 
      writeConstituent(&(lpa[ix]), std::cout);
  }
  */
  // Test out new utilities Eric requested
  //   getConfigsForAlias
  std::set<unsigned> keySet;

  bool ok = q.getConfigsForAlias("CAL_FLE", "100MeV", keySet);
  std::cout << "getConfigsForAlias(CAL_FLE, 100MeV..) found " <<
    keySet.size() << " configs" << std::endl << std::endl;

  keySet.clear();
  //   getConfigsForAncillary
  unsigned ancKey=29;
  ok = q.getConfigsForAncillary(ancKey, keySet);
  std::cout << "getConfigsForAncillary(26..) found " <<
    keySet.size() << " configs" << std::endl << std::endl;

  keySet.clear();
  //   getConfigsForVote
  voteKey=100;
  ok = q.getConfigsForVote(voteKey, keySet);
  std::cout << "getConfigsForVote(100..) found " <<
    keySet.size() << " configs" << std::endl << std::endl;
  

  //   getVotesForPrecinct
  std::vector<MOOT::VoteInfo> vInfo;
  ok = q.getVotesForPrecinct("CAL_Mode", vInfo);
  std::cout << "getVotesForPrecint(CAL_Mode..) found " << vInfo.size()
            << " votes " << std::endl;


  // getConstituentsBySbs
  std::vector<MOOT::ConstitInfo> conInfo;
  unsigned sbsExportKey = 6241;
  unsigned nAll = q.getConstituentInfoBySbs(sbsExportKey, conInfo);
  std::cout << "getConstituentsBySbs for sbs_key=" << sbsExportKey
            << " returns " << nAll << " constituents." << std::endl;
   
  unsigned schemaId = 225;  // gamma

  conInfo.clear();

  unsigned nGamma = q.getConstituentInfoBySbs(sbsExportKey, conInfo, schemaId);
  std::cout << "getConstituentsBySbs for sbs_key=" << sbsExportKey
            << " schema id =" << schemaId << " returns " << nGamma
            << " constituents." << std::endl;

  // Cache contents of  filter assoc file
  //  std::string assocFile("/u/ey/jrb/junk/filterAssoc.xml");
  //  q.parseFilterParm(assocFile);
//  unsigned cfgKey=2595;
  unsigned cfgKey=2587;
  //  q.cacheFilterConfig(cfgKey);
  conInfo.clear();
  unsigned nFound = q.getActiveFilters(cfgKey, conInfo);
  std::cout << "For cfgKey " << cfgKey << " found " << nFound << " filters"
            << std::endl;
  if (nFound) {
    writeConstituent(&conInfo[0], std::cout);
  }

  conInfo.clear();
  nFound = q.getActiveFilters(cfgKey, conInfo, MOOT::LPA_MODE_NORMAL);
  std::cout << "For cfgKey " << cfgKey << " and mode = NORMAL found " 
            << nFound << " filters" << std::endl;
  if (nFound) {
    writeConstituent(&conInfo[0], std::cout);
  }
  conInfo.clear();
  nFound = q.getActiveFilters(cfgKey, conInfo, MOOT::LPA_MODE_GRB1);
  std::cout << "For cfgKey " << cfgKey << " and mode = GRB1 found " 
            << nFound << " filters" << std::endl;
  if (nFound) {
    writeConstituent(&conInfo[0], std::cout);
  }

  std::string handlerName;
  for (unsigned hId = 0; hId < 6; hId++) {
    const MOOT::ConstitInfo* pInfo = 
      q.getActiveFilter(cfgKey, MOOT::LPA_MODE_NORMAL, hId, handlerName);
    if (pInfo) {
      std::cout << "For mode NORMAL, handler id " << hId 
                << " handler name is " << handlerName << std::endl;
      writeConstituent(pInfo, std::cout);
    }
  }

  MOOT::AcqSummaryInfo* pAcqInfo = q.getAcqSummaryInfoByKey(7);
  if (pAcqInfo) {
    std::cout << "Result of getAcqSummaryInfoByKey(7):" << std::endl;
    writeAcqInfo(pAcqInfo, std::cout);
  }
  else std::cout << "No info for key = 7" << std::endl << std::endl;

  unsigned startedAt = 223083791;
  pAcqInfo = q.getAcqSummaryInfo(startedAt);
  if (pAcqInfo) {
    std::cout << "Result of getAcqSummaryInfo(223083791):" << std::endl;
    writeAcqInfo(pAcqInfo, std::cout);
  }
  else std::cout << "No info for started_at = 223083791" 
                 << std::endl << std::endl;

  facilities::Timestamp t("2008-1-27");
  pAcqInfo = q.getAcqSummaryInfo(t);
  if (pAcqInfo) {
    std::cout << "Result of getAcqSummaryInfo(time) where time is 2008-1-27:" 
              << std::endl;
    writeAcqInfo(pAcqInfo, std::cout);
  }
  else std::cout << "No info for time 2008-1-27" << std::endl << std::endl;

  unsigned configKey=138;
  unsigned masterKey = q.getMasterKey(configKey);
  std::cout << "For Config #" << configKey << " Retrieved master key "
            << masterKey << std::endl;

  return 0;
}
void writeInfo(MOOT::ConfigInfo* pInfo, std::ostream& out) {
  out << "Key = " << pInfo->getKey() << std::endl;
  out << "Name = " << pInfo->getName() << std::endl;
  out << "Alg. = " << pInfo->getAlg() << std::endl;
  out << "Step = " << pInfo->getStep() << std::endl;
  out << "Description = " << pInfo->getDescrip() << std::endl;
  out << "Status = " << pInfo->getStatus() << std::endl;
  out << "active_state = " << pInfo->getActive() << std::endl;
  out << "Mode = " << pInfo->getMode() << std::endl;
  out << "Creation time = " << pInfo->getCreationTime() << std::endl;
  out << "Vote key = '" << pInfo->getVoteKey() << "'" << std::endl;
  out << "sbs export key = '" << pInfo->getSbsKey() << "'" << std::endl
            << std::endl;
}

void writeParmOff(const std::vector<MOOT::ParmOffline>& parm) {
  for (unsigned ix = 0; ix < parm.size(); ix++) {
    std::cout << "Parameter key = " << parm[ix].getKey() << std::endl;
    std::cout << "Class key = " << parm[ix].getClassFk() 
              <<  " with class name " << parm[ix].getClass()  << std::endl;
    std::cout << "Relative file path: " << parm[ix].getSrc()  << std::endl;
    std::cout << "Source format: " << parm[ix].getSrcFmt() 
              << "   Status: " << parm[ix].getStatus()  
              << std::endl << std::endl;
  }

}

void writeConstituent(const MOOT::ConstitInfo* pInfo, std::ostream& out) {
  out << "Constituent info: " << std::endl;
  out << "prim_key = " << pInfo->getKey() << std::endl;
  out << "name = " << pInfo->getName() << std::endl;
  out << "package = " << pInfo->getPkg() << std::endl;
  out << "package version = " << pInfo->getPkgVersion() << std::endl;
  out << "dir = " << pInfo->getDir() << std::endl;
  out << "fmx_path = " << pInfo->getFmxPath() << std::endl;
  out << "src_path = " << pInfo->getSrcPath() << std::endl;
  out << "fsw id = " << pInfo->getFswId() << std::endl;
  out << "status = " << pInfo->getStatus() << std::endl;
  if  ((pInfo->getSchemaId()).size() > 0 ) {
    out << "schema id = " << pInfo->getSchemaId() << std::endl;
    out << "schema version id = " << pInfo->getSchemaVersionId() 
        << std::endl;
    out << "instance id = " << pInfo->getInstanceId() << std::endl;
  }
  out << std::endl;
}

void writeAcqInfo(const MOOT::AcqSummaryInfo* pInfo, std::ostream& out) {
  out << "Acq summary info: " << std::endl;
  out << "acq_summary_key = " << pInfo->getKey() << std::endl;
  out << "scid = " << pInfo->getScid() << std::endl;
  out << "started_at = " << pInfo->getStartedAt() << std::endl;
  out << "config_fk = " << pInfo->getConfigKey() << std::endl;
  out << "type = " << pInfo->getAcqType() << std::endl;
  out << "analaysis = " << pInfo->getAnalysis() << std::endl;
  out << "status = " << pInfo->getStatus() << std::endl;
  out << "nevts = "  << pInfo->getNevts() << std::endl;
  out << "evtutc0 = "  << pInfo->getEvtutc0() << std::endl;
  out << "evtutc1 = "  << pInfo->getEvtutc1() << std::endl;
  out << "hwKey = "  << pInfo->getHwKey() << std::endl;
  out << "swKey = "  << pInfo->getSwKey() << std::endl;
  out << "creator = "  << pInfo->getCreator() << std::endl;
  out << "creation_time = "  << pInfo->getCreationTime() << std::endl;
  out << "update_time = "  << pInfo->getUpdateTime() << std::endl;
  out << std::endl;
}
