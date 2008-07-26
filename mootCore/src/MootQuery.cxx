//  $Header$

#include <string>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <set>
#include "mootCore/MootQuery.h"
#include "mootCore/MoodConnection.h"
#include "mootCore/DbUtil.h"
#include "rdbModel/Rdb.h"
#include "rdbModel/Tables/Table.h"
#include "rdbModel/Db/Connection.h"
#include "rdbModel/Db/ResultHandle.h"
#include "facilities/Util.h"
#include "facilities/Timestamp.h"
#include "mootCore/DbUtilException.h"

namespace {
  unsigned getAncsFromP(const std::vector<unsigned>& parms,
                        rdbModel::Rdb* rdb,
                        std::vector<unsigned>& ancs) {
    unsigned cnt;
    for (unsigned parmIx = 0; parmIx < parms.size(); parmIx++) {
      std::string parmKeyStr;
      facilities::Util::utoa(parms[parmIx], parmKeyStr);
      std::string where(" WHERE Parameter_fk='");
      where += parmKeyStr + std::string("'");
      // Get all anc keys for this parameter in ascending order
      cnt += MOOT::DbUtil::getKeys(ancs, rdb, "Parameters_to_Ancillary", 
                                   "Ancillary_fk", where, 0, true);
    }
    return cnt;
  }

}

namespace MOOT {
  MootQuery::MootQuery(MoodConnection* mood) : m_rdb(0), m_mood(mood), 
                                               m_ownConnection(false),
                                               m_archive("")  {
    m_parmClasses.clear();
    if (!mood) {
      m_mood = new MoodConnection();
      if (!m_mood) {
        throw std::runtime_error("MootQuery constructor: can't construct MoodConnection");
      }
      else m_ownConnection = true;
    }
    m_rdb = m_mood->getRdb();
    if (!m_rdb) {
      m_mood = 0;
      throw 
        std::runtime_error("MootQuery constructor: can't get rdbModel::Rdb ");
    }
    m_archive = std::string("$(MOOT_ARCHIVE)/");
    m_dbg = m_mood->isDebug();

    try {
      int nVar = facilities::Util::expandEnvVar(&m_archive);
      if (nVar == 1 ) return;
    }
    catch (facilities::Untranslatable ex) {
       { // we're toast
         std::cerr << "Missing definition for env. variable MOOT_ARCHIVE " 
                   << std::endl;
         std::cerr << "and can't do much without it.  Bye." << std::endl;
         std::cerr.flush();
         throw ex;
       }
    } 

  }
  MootQuery::~MootQuery() {
    if (m_ownConnection) delete m_mood;
  }

  int MootQuery::classKey(const std::string& table, const std::string& name) {
    std::string getCol = m_rdb->getTable(table)->getPrimaryKeyCol();
    std::string keyStr;
    try {
      keyStr = DbUtil::getColumnValue(m_rdb, table, getCol, "name", name);
    }
    catch (std::exception ex) {
      std::cerr << ex.what() << std::endl;
      std::cerr.flush();
      return -1;
    }
    if (keyStr.size() == 0) return -1;
    return facilities::Util::stringToInt(keyStr);
  }

  std::string MootQuery::classStr(const std::string& table, unsigned key) {
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    std::string keyCol = m_rdb->getTable(table)->getPrimaryKeyCol();
    try {
      return DbUtil::getColumnValue(m_rdb, table, "name", keyCol, keyStr);
    }
    catch (std::exception ex) {
      std::cerr << ex.what() << std::endl;
      std::cerr.flush();
      throw ex;
    }
  }

  AcqSummaryInfo* MootQuery::getAcqSummaryInfoByKey(unsigned key) {
    std::string where(" WHERE acq_summary_key='");
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    where += keyStr + std::string("'");
    return getAcqSummaryWhere(where);
  }

  AcqSummaryInfo* MootQuery::getAcqSummaryInfo(unsigned startedAt,
                                               unsigned scid) {
    std::string scidStr;
    facilities::Util::utoa(scid, scidStr);
    std::string startedStr;
    facilities::Util::utoa(startedAt, startedStr);
    std::string where(" WHERE scid='");
    where += scidStr + std::string("' AND started_at='");
    where += startedStr + std::string("'");
    return getAcqSummaryWhere(where);
  }

  AcqSummaryInfo* 
  MootQuery::getAcqSummaryInfo(const facilities::Timestamp& eventTime,
                               unsigned scid) {
    std::string scidStr;
    facilities::Util::utoa(scid, scidStr);
   
    std::string timStr = eventTime.getString();
    std::string where(" WHERE scid='");
    where += scidStr + std::string("' AND evtutc0 <='") + timStr;
    where += std::string("' AND evtutc1 >='") + timStr + std::string("'");
    return getAcqSummaryWhere(where);
  }

  AcqSummaryInfo* MootQuery::getAcqSummaryWhere(const std::string& where) {
    std::vector<std::string> getCols;
    std::vector<std::string> noCols;
    getCols.reserve(15);

    getCols.push_back("acq_summary_key");
    getCols.push_back("scid");
    getCols.push_back("started_at");
    getCols.push_back("config_fk");
    getCols.push_back("type");
    getCols.push_back("analysis");
    getCols.push_back("status");
    getCols.push_back("nevts");
    getCols.push_back("evtutc0");
    getCols.push_back("evtutc1");
    getCols.push_back("hw_key");
    getCols.push_back("sw_key");
    getCols.push_back("creator");
    getCols.push_back("creation_time");
    getCols.push_back("update_time");

    rdbModel::ResultHandle* res = 0;

    try {
      res = m_rdb->getConnection()->select("Acq_summary", 
                                           getCols, noCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getAcqSummaryInfo, where = " << where
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }
    int n = res->getNRows();
    if (!n) {
      delete res;
      return 0;
    }
    
    std::vector<std::string>  fields;
    res->getRow(fields, 0);
    delete res;
    return new AcqSummaryInfo(fields[0], fields[1], fields[2], fields[3],
                              fields[4], fields[5], fields[6], fields[7],
                              fields[8], fields[9], fields[10], fields[11],
                              fields[12], fields[13], fields[14]);
  }

  AncAliasInfo* MootQuery::getAncAliasInfo(unsigned key) {
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    std::string where(" WHERE ancillary_aliases_key='");
    where += keyStr + std::string("'");

    std::vector<std::string> getCols;
    std::vector<std::string> noCols;
    getCols.reserve(5);
    getCols.clear(); noCols.clear();

    getCols.push_back("name");
    getCols.push_back("aclass_fk");
    getCols.push_back("ancillary_fk");
    getCols.push_back("last_modified_time");
    getCols.push_back("last_modifier");

    rdbModel::ResultHandle* res = 0;

    try {
      res = m_rdb->getConnection()->select("Ancillary_aliases", 
                                           getCols, noCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getAncAliasInfo, key = " << key
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }
    int n = res->getNRows();
    if (!n) {
      delete res;
      return 0;
    }
    
    std::vector<std::string>  fields;
    res->getRow(fields, 0);
    delete res;

    // translate anc class key to anc class name
    std::string ancClassName =
      DbUtil::getColumnValue(m_rdb, "Ancillary_class", "name", 
                             "Ancillary_class_key", fields[1]);
    return new AncAliasInfo(keyStr, fields[0], ancClassName, fields[1],
                            fields[2],fields[3],fields[4]);
  }


  unsigned MootQuery::getAncClasses(std::vector<std::string>& names) {
    std::string where("");
    names.clear();
    int n = DbUtil::getAllWhere(m_rdb, "Ancillary_class", "name", where,
                                names);
    return (n >= 0) ? 1 : 0;
  }

  AncInfo* MootQuery::getAncInfo(unsigned key) {
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    std::string where(" WHERE Ancillary_key='");
    where += keyStr + std::string("'");

    std::vector<std::string> getCols;
    std::vector<std::string> noCols;
    getCols.reserve(10);
    getCols.clear(); noCols.clear();

    getCols.push_back("class_fk");
    getCols.push_back("instrument");
    getCols.push_back("source");
    getCols.push_back("source_fmt");
    getCols.push_back("GMT_create_time");
    getCols.push_back("creation_time");
    getCols.push_back("creator");
    getCols.push_back("quality");
    getCols.push_back("description");
    getCols.push_back("status");

    rdbModel::ResultHandle* res = 0;

    try {
      res = m_rdb->getConnection()->select("Ancillary", 
                                           getCols, noCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getAncInfo, key = " << key
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }
    int n = res->getNRows();
    if (!n) {
      delete res;
      return 0;
    }
    
    std::vector<std::string>  fields;
    res->getRow(fields, 0);
    delete res;
    // translate anc class key to anc class name
    std::string ancClassName =
      DbUtil::getColumnValue(m_rdb, "Ancillary_class", "name", 
                             "Ancillary_class_key", fields[0]);
    return new AncInfo(keyStr, ancClassName, fields[0], fields[1],
                       fields[2],fields[3],fields[4], fields[5],
                       fields[6],fields[7], fields[8],fields[9]);
  }

  unsigned MootQuery::getAncsFromParms(const std::vector<unsigned>& parms,
                                       std::vector<unsigned>& ancs) {
    return getAncsFromP(parms, m_rdb, ancs);
  }

  bool MootQuery::getConfigAncsRequest(unsigned configKey,
                                       std::vector<unsigned>& ancKeys) {
    std::vector<unsigned> parmKeys;
    bool ok = getConfigParmsRequest(configKey, parmKeys);
    if (!ok) return ok;
    getAncsFromP(parmKeys, m_rdb, ancKeys);
    return true;
  }

  bool MootQuery::getConfigAncsUsed(unsigned configKey,
                                    std::vector<unsigned>& ancKeys) {
    std::vector<unsigned> parmKeys;
    bool ok = getConfigParmsUsed(configKey, parmKeys);
    if (!ok) return ok;
    getAncsFromP(parmKeys, m_rdb, ancKeys);
    return true;
  }

  bool MootQuery::getConfigFmxPaths(unsigned configKey,
                                    std::vector<std::string>& paths) {
    std::string keyStr;
    facilities::Util::utoa(configKey, keyStr);
    return getConfigFmxPaths(keyStr, paths);
  }

  bool MootQuery::getConfigFmxPaths(const std::string& configKey,
                                    std::vector<std::string>& paths) {
    std::vector<std::string> inputKeys;
    bool ret = true;

    if (!getConfigInputs(configKey, inputKeys)) return false;

    for (unsigned i = 0; i < inputKeys.size(); i++ ) {
      std::string status;
      try {
        status = DbUtil::getColumnValue(m_rdb, "FSW_inputs", "status",
                                        "FSW_input_key", inputKeys[i]);
      }
      catch (std::exception ex) {
        std::cerr << ex.what() << std::endl;
        std::cerr.flush();
        return false;
      }
      if (status == std::string("added")) {
        try {
          std::string path = 
            DbUtil::getColumnValue(m_rdb, "FSW_inputs", "fmx_path", 
                                   "FSW_input_key", inputKeys[i]);
          paths.push_back(path);
        }
        catch (std::exception ex) {
          std::cerr << ex.what() << std::endl;
          std::cerr.flush();
          return false;
        }
      }

      else {
        std::cerr << "Config " << configKey 
                  << "includes unusable FSW input with MOOT id " 
                  << inputKeys[i]  << std::endl;
        std::cerr.flush();
        ret = false;
      }
    }
    return ret;
  }

  bool MootQuery::getConfigFmxPathsByAlg(const std::string& algName,
                                         const std::string& algStep,
                                         std::vector<std::string>& paths) {
    std::string where(" WHERE algorithm='");
    where += algName + std::string("' AND alg_step='") + algStep +
      std::string("' AND status='CREATED' ORDER BY config_key DESC" );
    std::string key;
    try {
      key = DbUtil::getColumnWhere(m_rdb, "Configs", "config_key", where, 
                                   false);
    }
    catch (std::exception ex) {
      std::cerr << ex.what() << std::endl;
      std::cerr.flush();
      return false;
    }
    if (key.size() == 0 ) return false;
    return getConfigFmxPaths(key, paths);
  }

  bool MootQuery::getConfigFmxPathsByAlg(const std::string& algName,
                                         unsigned algStep,
                                         std::vector<std::string>& paths)
  {
    std::string stepStr;
    facilities::Util::utoa(algStep, stepStr);
    return getConfigFmxPathsByAlg(algName, stepStr, paths);
  }

  bool MootQuery::getConfigFmxPathsByName(const std::string& configName,
                                          std::vector<std::string>& paths) {
    std::string where(" WHERE name='");
    where += configName +
      std::string("' AND status='CREATED' ORDER BY config_key DESC" );
    std::string key;
    try {
      key = DbUtil::getColumnWhere(m_rdb, "Configs", "config_key",
                                   where, false);
    }
    catch (std::exception ex) {
      std::cerr << ex.what() << std::endl;
      std::cerr.flush();
      return false;
    }
    if (key.size() == 0 ) return false;
    return getConfigFmxPaths(key, paths);
  }


  bool MootQuery::getConfigFswInfo(unsigned configKey,
                                   std::vector<FswInfo>& info, bool clear) {
    std::string keyStr;
    facilities::Util::utoa(configKey, keyStr);
    return getConfigFswInfo(keyStr, info, clear);
  }

  bool MootQuery::getConfigFswInfo(const std::string& configKey,
                                   std::vector<FswInfo>& info, bool clear) {

    std::vector<std::string> inputKeys;
    if (clear) info.clear();

    if (!getConfigInputs(configKey, inputKeys)) return false;
    rdbModel::StringVector getCols;
    rdbModel::StringVector noCols;
    getCols.reserve(6);
    getCols.push_back("fmx_path");
    getCols.push_back("source");
    getCols.push_back("description");
    getCols.push_back("FSW_id");
    getCols.push_back("class_fk");
    getCols.push_back("fmx_db");

    for (unsigned i = 0; i < inputKeys.size(); i++ ) {
      std::string where(" WHERE FSW_input_key=");
      where += inputKeys[i];
      rdbModel::ResultHandle* res = 0;
      try {
        res = m_rdb->getConnection()->select("FSW_inputs", getCols, noCols,
                                             where);
      }
      catch (std::exception ex) {
        std::cerr << "MootQuery::getConfigFswInfo, inputKey= " << inputKeys[i]
                  << " SQL error: "  << ex.what() << std::endl;
        std::cerr.flush();
        if (res) delete res;
        return false;
      }
      std::vector<std::string>selFields;
      res->getRow(selFields);
      unsigned id = facilities::Util::stringToUnsigned(selFields[3]);
      unsigned classFk = facilities::Util::stringToUnsigned(selFields[4]);
      FswInfo f(selFields[0], selFields[1], selFields[2], selFields[5],
                id, classFk);
      info.push_back(f);
      delete res;
    }
    return true;
  }

  bool MootQuery::getConfigFswInfoByAlg(const std::string& algName,
                                        unsigned algStep,
                                        std::vector<FswInfo>& info,
                                        bool clear) {
    std::string stepStr;
    facilities::Util::utoa(algStep, stepStr);
    return getConfigFswInfoByAlg(algName, stepStr, info, clear);
  }

  bool MootQuery::getConfigFswInfoByAlg(const std::string& algName,
                                        const std::string& algStep,
                                        std::vector<FswInfo>& info,
                                        bool clear) {
    std::string where(" WHERE algorithm='");
    where += algName + std::string("' AND alg_step='") + algStep +
      std::string("' AND status='CREATED' ORDER BY config_key DESC" );
    std::string key;
    try {
      key = DbUtil::getColumnWhere(m_rdb, "Configs", "config_key",
                                   where, false);
    }
    catch (std::exception ex) {
      std::cerr << ex.what() << std::endl;
      std::cerr.flush();
      return false;
    }

    if (key.size() == 0 ) return false;
    return getConfigFswInfo(key, info, clear);
  }

  bool MootQuery::getConfigFswInfoByName(const std::string& configName,
                                         std::vector<FswInfo>& info, 
                                         bool clear) {
    std::string where(" WHERE name='");
    where += configName +
      std::string("' AND status='CREATED' ORDER BY config_key DESC" );
    std::string key;
    try {
      key = DbUtil::getColumnWhere(m_rdb, "Configs", "config_key",
                                   where, false);
    }
    catch (std::exception ex) {
      std::cerr << ex.what() << std::endl;
      std::cerr.flush();
      return false;
    }
    if (key.size() == 0 ) return false;
    return getConfigFswInfo(key, info, clear);
  }

  ConfigInfo* MootQuery::getConfigInfo(unsigned key) {
    rdbModel::StringVector getCols;
    rdbModel::StringVector orderCols;
    orderCols.clear();
    std::string where(" WHERE config_key= '");
    std::string keystr;
    facilities::Util::utoa(key, keystr);
    where += keystr + "'";
    
    getCols.reserve(12);
    getCols.push_back("config_key");
    getCols.push_back("name");
    getCols.push_back("algorithm");
    getCols.push_back("alg_step");
    getCols.push_back("description");
    getCols.push_back("status");
    getCols.push_back("active_state");
    getCols.push_back("verify_state");
    getCols.push_back("mode");
    getCols.push_back("creation_request_time");
    getCols.push_back("vote_fk");
    getCols.push_back("sbs_key");
    rdbModel::ResultHandle* res = 0;
    try {
      res = m_rdb->getConnection()->select("Configs", getCols, orderCols,
                                           where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getConfigInfo "
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      return 0;
    }
    int n = res->getNRows();
    if (!n) {
      delete res;
      return 0;
    }
    
    std::vector<std::string>selFields;
    res->getRow(selFields, 0);
    ConfigInfo* c = 
      new ConfigInfo(selFields[0], selFields[1], selFields[2], selFields[3],
                     selFields[4], selFields[5], selFields[6], selFields[7],
                     selFields[8], selFields[9], selFields[10], selFields[11]);
    delete res;
    return c;
  }

  unsigned MootQuery::getConfigInfo(std::vector<ConfigInfo>& info,
                                    const std::string& status,
                                    const std::string& instr,
                                    const std::string& mode) {
    std::string where("");
    if (status.size()) {
      where = std::string(" WHERE status='") + status + std::string("' ");
    }
    if (instr.size() ) {
      if (where.size()) {
        where += std::string(" AND instrument='") + instr + std::string("' ");
      }
      else {
        where = std::string(" WHERE instrument='") + instr + std::string("' ");
      }
    }

    if (mode.size() ) {
      if (where.size()) {
        where += std::string(" AND mode='") + mode + std::string("' ");
      }
      else {
        where = std::string(" WHERE mode='") + mode + std::string("' ");
      }
    }
    rdbModel::StringVector getCols;
    rdbModel::StringVector orderCols;

    orderCols.push_back("config_key");
    getCols.reserve(11);
    getCols.push_back("config_key");
    getCols.push_back("name");
    getCols.push_back("algorithm");
    getCols.push_back("alg_step");
    getCols.push_back("description");
    getCols.push_back("status");
    getCols.push_back("active_state");
    getCols.push_back("verify_state");
    getCols.push_back("mode");
    getCols.push_back("creation_request_time");
    getCols.push_back("vote_fk");
    getCols.push_back("sbs_key");

    rdbModel::ResultHandle* res = 0;
    try {
      res = m_rdb->getConnection()->select("Configs", getCols, orderCols,
                                           where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getConfigInfo "
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      return 0;
    }
    int n = res->getNRows();

    info.reserve(info.size() + n);
    std::vector<std::string>selFields;
    
    for (int i = 0; i < n; i++) {
      res->getRow(selFields, i);
      ConfigInfo info1(selFields[0], selFields[1], selFields[2], selFields[3],
                       selFields[4], selFields[5], selFields[6], selFields[7],
                       selFields[8], selFields[9], selFields[10], selFields[11]);
      info.push_back(info1);
    }
    delete res;
    return n;
  }


  bool MootQuery::getConfigInputs(const std::string& configKeyStr,
                                  std::vector<std::string>& inputKeys,
                                  bool goodStatus) {
    if (goodStatus) {
      std::string configStatus;
      try {
        configStatus = 
          DbUtil::getColumnValue(m_rdb, "Configs", "status", "config_key",
                                 configKeyStr);
      }
      catch (std::exception ex) {
        std::cerr << ex.what() << std::endl;
        std::cerr.flush();
        return false;
      }
      if (configStatus != std::string("CREATED") ) {
        std::cerr << "MootQuery::getConfigInputs:  Config "  
                  << configKeyStr << " does not have good (CREATED) status " 
                  << std::endl;
        std::cerr.flush();
        return false;
      }
    }
    inputKeys.clear();
    std::string where(" WHERE config_fk='");
    where += configKeyStr + "' ";
    // int nFetched = 
    DbUtil::getAllWhere(m_rdb, "Configs_to_FSW", "FSW_fk", 
                                       where, inputKeys);
    return true;
  }

  bool MootQuery::getConfigInputs(unsigned configKey, 
                                  std::vector<std::string>& inputKeys,
                                  bool goodStatus) {

    std::string keyStr;
    facilities::Util::utoa(configKey, keyStr);
    return getConfigInputs(keyStr, inputKeys, goodStatus);
  }

  unsigned MootQuery::getConfigKeysByAlg(std::vector<unsigned>& keys, 
                                         const std::string& alg,
                                         unsigned step,
                                         const std::string& status,
                                         const std::string& instr) {
    std::string where(" WHERE algorithm='");
    where += alg + std::string("'");
    if (step > 0) {
      std::string stepStr;
      facilities::Util::utoa(step, stepStr);
      where += std::string(" AND alg_step='")+ stepStr + "'";
    }
    if (status.size()) {
      where += std::string(" AND status='" ) + status +
        std::string("'");
    }
    if (instr.size()) {
      where  += std::string(" AND instrument='" ) + instr +
        std::string("'");
    }
    return
      DbUtil::getKeys(keys, m_rdb, "Configs", "config_key", where,
                      0, true);   // keys in ascending order

  }

  unsigned MootQuery::getConfigKeysByName(std::vector<unsigned>& keys, 
                                          const std::string& name, 
                                          const std::string& status,
                                          const std::string& instr) {
    std::string where(" WHERE name='");
    where += name + std::string("'");
    if (status.size()) {
      where += std::string(" AND status='" ) + status +
        std::string("'");
    }
    if (instr.size()) {
      where  += std::string(" AND instrument='" ) + instr +
        std::string("'");
    }

    // Ask for keys in ascending order
    return DbUtil::getKeys(keys, m_rdb, "Configs", "config_key", where,
                           0, true);
  }

  bool MootQuery::getConfigParmsRequest(unsigned configKey, 
                                        std::vector<unsigned>& parameterKeys){
    std::string where("WHERE config_fk ='");
    std::string configKeyStr;
    facilities::Util::utoa(configKey, configKeyStr);
    where += configKeyStr + std::string("'");
    DbUtil::getKeys(parameterKeys, m_rdb, "Configs_to_Parameters",
                    "Parameter_fk", where, 0, true);
    return true;
  }

  bool MootQuery::getConfigParmsUsed(unsigned configKey, 
                                     std::vector<unsigned>& parameterKeys) {
    std::vector<std::string> inputKeys;
    bool ok = getConfigInputs(configKey, inputKeys);
    if (!ok) return false;
    if (!inputKeys.size()) return ok;

    // For each input..
    for (unsigned i = 0; i < inputKeys.size(); i++) {
      std::string where(" WHERE FSW_fk='");
      where += inputKeys[i] + "'";
      // get all parameter keys assoc. with it in ascending order 
      DbUtil::getKeys(parameterKeys, m_rdb, 
                                      "FSW_to_Parameters",
                                      "Parameter_fk", where, 0,
                                      true);
    }
    std::sort(parameterKeys.begin(), parameterKeys.end());
    return true;
  }
  bool MootQuery::getConfigsForAlias(const std::string& precinct, 
                                     const std::string& alias,
                                     std::set<unsigned>& keys) {
    // get key corresponding to precinct name
    std::string where(" WHERE name='");
    where += precinct + std::string("'");
    std::string prKeyStr = 
      DbUtil::getColumnWhere(m_rdb, "Precincts", "precinct_key", where,
                             false);
    if (prKeyStr.empty()) return false;
    where = std::string(" WHERE Container_Precinct.v_alias='");
    where += alias + std::string("' AND Container_Precinct.precinct_fk ='");
    where += prKeyStr + 
      std::string("' AND Configs.vote_fk=Container_Precinct.ctn_fk");
    std::vector<std::string> getCols;
    getCols.push_back("Configs.config_key");
    rdbModel::ResultHandle* res = 
      m_rdb->getConnection()->select("Container_Precinct,Configs",
                                     getCols, getCols, where);

    if (!res) return false;  // shouldn't happen even if none are found
    int nRows = res->getNRows();
    std::vector<std::string> vals;
    std::set<unsigned>::iterator hint=keys.begin();
    for (unsigned ix = 0; ix < nRows; ix++) {
      res->getRow(vals, ix);
      unsigned key = facilities::Util::stringToUnsigned(vals[0]);
      hint = keys.insert(hint, key);
    }
    delete res;
    return true;
  }

  bool MootQuery::getConfigsForAncillary(unsigned ancKey, 
                                         std::set<unsigned>& keys) {
    using facilities::Util;
    /*
         select config_fk from 
         Configs_to_Parameters,Parameters_to_Ancillary where 
         Ancillary_fk=29 and Configs_to_Parameters.Parameter_fk = 
         Parameters_to_Ancillary.Parameter_fk;
    */
    std::string ancKeyStr;
    std::string where(" WHERE Ancillary_fk='");
    Util::utoa(ancKey, ancKeyStr);
    where += ancKeyStr + 
      std::string("' AND Configs_to_Parameters.Parameter_fk=Parameters_to_Ancillary.Parameter_fk ");

    std::vector<std::string> getCols;
    getCols.push_back("config_fk");
    rdbModel::ResultHandle* res = 
      m_rdb->getConnection()->select("Configs_to_Parameters,Parameters_to_Ancillary",
                                     getCols, getCols, where);

    if (!res) return false;  // shouldn't happen even if none are found
    int nRows = res->getNRows();
    std::vector<std::string> vals;
    std::set<unsigned>::iterator hint=keys.begin();
    for (unsigned ix = 0; ix < nRows; ix++) {
      res->getRow(vals, ix);
      unsigned key = Util::stringToUnsigned(vals[0]);
      hint = keys.insert(hint, key);
    }
    delete res;
    return true;
  }

  bool MootQuery::getConfigsForVote(unsigned voteKey,
                                    std::set<unsigned>& keys) {
    using facilities::Util;
    /*
      select config_fk from  Parameters,Configs_to_Parameters where 
       vote_fk=120 and parm_key=Parameter_fk;
    */
    std::string voteKeyStr;
    Util::utoa(voteKey, voteKeyStr);
    std::string where(" WHERE vote_fk='");
    where += voteKeyStr + 
      std::string("' AND Parameters.parm_key=Configs_to_Parameters.Parameter_fk");
    std::vector<std::string> getCols;
    getCols.push_back("config_fk");
    rdbModel::ResultHandle* res = 
      m_rdb->getConnection()->select("Parameters,Configs_to_Parameters",
                                     getCols, getCols, where);

    if (!res) return false;  // shouldn't happen even if none are found
    int nRows = res->getNRows();
    std::vector<std::string> vals;
    std::set<unsigned>::iterator hint=keys.begin();
    for (unsigned ix = 0; ix < nRows; ix++) {
      res->getRow(vals, ix);
      unsigned key = Util::stringToUnsigned(vals[0]);
      hint = keys.insert(hint, key);
    }
    delete res;
    return true;
  }

  ConstitInfo* MootQuery::getConstituentByFswId(unsigned fswId) {
    std::string fswIdStr;
    facilities::Util::utoa(fswId, fswIdStr);
    std::string where(" WHERE FSW_id='");
    where += fswIdStr + std::string("'");
    return getConstituentWhere(where);
  }

  ConstitInfo* MootQuery::getConstituentInfo(unsigned key) {
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    std::string where(" WHERE prim_key= '");
    where += keyStr + std::string("'");
    return getConstituentWhere(where);
  }

  unsigned 
  MootQuery::getConstituentInfoBySbs(unsigned sbsFmxExportKey, 
                                     std::vector<ConstitInfo>& info,
                                     int schemaId) {
    using facilities::Util;
    std::string where(" WHERE FSW_id ='");
    std::string keyStr;
    Util::utoa(sbsFmxExportKey, keyStr);
    where += keyStr + std::string("'");

    std::string sbsFk = DbUtil::getColumnWhere(m_rdb, "FSW_inputs", 
                                               "FSW_input_key", where,
                                               false);
    if (sbsFk.empty()) return 0;

    /// Constraints are
    //    SBS_to_Constituents.sbs_fk = sbsFk (computed above)
    //    Constituents.prim_key = SBS_to_Constituents.Constituents_fk
    //    Constituents.schema_id = supplied value, if not -1
    std::vector<std::string> constitKeys;
    std::string tbls("SBS_to_Constituents");


    where = std::string(" WHERE sbs_fk='");
    where += sbsFk + std::string("'");
    if (schemaId >= 0) {
      tbls += std::string(",Constituents");
      std::string schemaIdStr;
      Util::utoa(schemaId, schemaIdStr);
      where +=
        std::string(" AND Constituents.prim_key=SBS_to_Constituents.Constituents_fk AND Constituents.schema_id='");
      where += schemaIdStr + std::string("'");
    }
    DbUtil::getAllWhere(m_rdb, tbls, "SBS_to_Constituents.Constituents_fk", 
                        where, constitKeys);

    info.reserve(constitKeys.size());
    for (unsigned ix = 0; ix < constitKeys.size(); ix++) {
      std::string where(" WHERE prim_key= '");
      where += constitKeys[ix] + std::string("'");
      ConstitInfo* pInfo =  getConstituentWhere(where);
      if (pInfo) info.push_back(*pInfo);
    }
    return info.size();
  }


  ConstitInfo* MootQuery::getConstituentWhere(const std::string& where) {

    std::vector<std::string> getCols;
    std::vector<std::string> noCols;
    getCols.reserve(12);
    getCols.clear(); noCols.clear();

    getCols.push_back("prim_key");
    getCols.push_back("name");
    getCols.push_back("pkg");
    getCols.push_back("ver");
    getCols.push_back("dir");
    getCols.push_back("fmx_path");
    getCols.push_back("src_path");
    getCols.push_back("FSW_id");
    getCols.push_back("status");
    getCols.push_back("schema_id");
    getCols.push_back("version_id");
    getCols.push_back("instance_id");
    rdbModel::ResultHandle* res = 0;

    try {
      res = m_rdb->getConnection()->select("Constituents", 
                                           getCols, noCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getConstituentWhere, where = " << where
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }
    int n = res->getNRows();
    if (!n) {
      delete res;
      return 0;
    }
    
    std::vector<std::string>  flds;
    res->getRow(flds, 0);
    delete res;
    return new ConstitInfo(flds[0], flds[1], flds[2], flds[3], flds[4],
                           flds[5], flds[6], flds[7], flds[8],
                           flds[9], flds[10], flds[11]);
  }

  unsigned MootQuery::getLastConfigKeyByAlg(const std::string& alg, 
                                            unsigned step,
                                            const std::string& status,
                                            const std::string& instr) {
    std::string where(" WHERE algorithm='");
    where += alg + std::string("'");
    if (step > 0) {
      std::string stepStr;
      facilities::Util::utoa(step, stepStr);
      where += std::string(" AND alg_step='")+ stepStr + "'";
    }
    if (status.size()) {
      where += std::string(" AND status='" ) + status +
        std::string("'");
    }
    if (instr.size()) {
      where  += std::string(" AND instrument='" ) + instr +
        std::string("'");
    }
    std::vector<unsigned> keys;
    unsigned nKeys = 
      DbUtil::getKeys(keys, m_rdb, "Configs", "config_key", where, 1);
    return (nKeys) ? keys[0] : 0;
  }

  unsigned MootQuery::getLastConfigKeyByName(const std::string& name, 
                                             const std::string& status,
                                             const std::string& instr) {
    std::string where(" WHERE name='");
    where += name + std::string("'");
    if (status.size()) {
      where += std::string(" AND status='" ) + status +
        std::string("'");
    }
    if (instr.size()) {
      where  += std::string(" AND instrument='" ) + instr +
        std::string("'");
    }
    std::vector<unsigned> keys;
    unsigned nKeys = 
      DbUtil::getKeys(keys, m_rdb, "Configs", "config_key", where, 1);
    return (nKeys) ? keys[0] : 0;
  }

  bool MootQuery::getLatcSrc(unsigned latcMasterKey,
                             std::vector<FileDescrip>& sources) {
    std::string keyStr;
    facilities::Util::utoa(latcMasterKey, keyStr);
    return getLatcSrc(keyStr, sources);
  }

  bool MootQuery::getLatcSrc(const std::string& latcMasterKey,
                             std::vector<FileDescrip>& sources) {
    std::string descrip;
    try {
      descrip = 
        DbUtil::getColumnValue(m_rdb, "FSW_inputs", "description",
                               "FSW_id", latcMasterKey);
    }
    catch (std::exception ex) {
      std::cerr << ex.what() << std::endl;
      std::cerr.flush();
      return false;
    }

    if (descrip.size() == 0) return false;

    // description for a latc master file looks like
    //   n_m_k_  ..and so forth, where n, m, k are (MOOT) keys for
    //   other rows in the table.
    unsigned pos = 0;
    while (pos   < descrip.size() ) {
      // parse out next key
      unsigned newPos = descrip.find('_', pos);
      std::string srcKey = std::string(descrip, pos, newPos - pos);
      //    get source and class_fk fields
      std::string path;
      try {
        path = DbUtil::getColumnValue(m_rdb, "FSW_inputs", "source", 
                                      "FSW_input_key", srcKey);
      }
      catch (std::exception ex) {
        std::cerr << ex.what() << std::endl;
        std::cerr.flush();
        return false;
      }

      if (path.size() == 0) {
        std::cerr << "MootQuery::getLatcSrc:  bad FSW_input key " 
                  << srcKey << std::endl;
      }
      else {
        path = m_archive + path;
        std::string where(" WHERE FSW_class_key=(SELECT class_fk from FSW_inputs  WHERE FSW_input_key =");
        where += srcKey + std::string(")");
        std::string srcType;
        try {
          srcType = DbUtil::getColumnWhere(m_rdb, "FSW_class", "name", where);
        }
        catch (std::exception ex) {
          std::cerr << ex.what() << std::endl;
          std::cerr.flush();
          return false;
        }
        sources.push_back(FileDescrip(path, srcType));
      }
      
      //    Prepend MOOT_ARCHIVE root to source to give abs. path
      //    look up name corresponding to class_fk in FSW_class table
      pos = newPos + 1;
    }
    return true;
  }

  // sbsKey arg.  refers to fmx logical key
  bool MootQuery::getLpaConstituents(unsigned sbsKey, 
                                     std::vector<ConstitInfo>& lpas) {
    return getPackageConstituents(sbsKey, "LPA_DB", lpas);
  }

  unsigned MootQuery::getMasterKey(unsigned configKey) {
    using facilities::Util;

    std::string configKeyStr;
    Util::utoa(configKey, configKeyStr);
    std::string masterClass = 
      DbUtil::getColumnValue(m_rdb, "FSW_class", "FSW_class_key", "name",
                             "latc_master");
    std::string where("WHERE FSW_inputs.class_fk='");
    where += masterClass + 
      std::string("' and FSW_input_key=FSW_fk AND config_fk='");
    where += configKeyStr + std::string("'");
    std::string masterKeyStr;
    try {
      masterKeyStr =
        DbUtil::getColumnWhere(m_rdb, "FSW_inputs,Configs_to_FSW", "FSW_id",
                               where);
    }
    catch (DbUtilException ex) {
      return 0;
    }
    return Util::stringToUnsigned(masterKeyStr);
  }

  bool MootQuery::getPackageConstituents(unsigned sbsKey, 
                                         const std::string& pkg,
                                         std::vector<ConstitInfo>& constits) {
    /*
      Find all constituents associated with
      the sbs key and with package = [pkg]

      select prim_key from SBS_to_Constituents where sbs_key = 
         (select FSW_input_key from FSW_inputs where FSW_id = [sbsKey] ) 
        and "LPA_DB" = (select pkg from Constituents where 
        Constituents.prim_key = SBS_to_Constituents.Constituents_fk)
    */
    std::string sbsKeyStr;
    facilities::Util::utoa(sbsKey, sbsKeyStr);
    std::vector<std::string> constitKeys;
    std::string 
      where("WHERE sbs_fk=(select FSW_input_key from FSW_inputs where FSW_id='");
    where += sbsKeyStr + 
      std::string("') and '") + pkg;
    where += std::string("'=(select pkg from Constituents where ");
    where += 
      std::string("Constituents.prim_key=SBS_to_Constituents.Constituents_fk)");
    int cnt = DbUtil::getAllWhere(m_rdb, "SBS_to_Constituents", "prim_key",
                                  where, constitKeys);
    if (cnt <= 0) return false;
    for (int ix=0; ix < cnt; ix++) {
      unsigned key = facilities::Util::stringToUnsigned(constitKeys[ix]);
      ConstitInfo* info = getConstituentInfo(key);
      constits.push_back(*info);
      delete info; 
    }
    return true;
  }

  unsigned MootQuery::getParameterClasses(std::vector<std::string>& names) {
    return getParmClasses(names);
  }
  unsigned MootQuery::getParmClasses(std::vector<std::string>& names,
                                     const std::string& precinct) {
    int nRet;
    if (precinct == std::string("*")) {
      nRet = DbUtil::getAllWhere(m_rdb, "Parameter_class", "name",
                                 "", names);
    }
    else {
      std::string where(" WHERE precinct_fk = (select precinct_key from Precincts where name='");
      where += precinct + std::string("')");
      nRet = DbUtil::getAllWhere(m_rdb, "Parameter_class", "name", where,
                                 names);
    }
    return nRet;
  }

  ParmInfo* MootQuery::getParmInfo(unsigned key) {
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    std::string where(" WHERE parm_key='");
    where += keyStr + std::string("'");

    std::vector<std::string> getCols;
    std::vector<std::string> noCols;
    getCols.reserve(10);
    getCols.clear(); noCols.clear();

    getCols.push_back("class_fk");
    getCols.push_back("vote_fk");
    getCols.push_back("instrument");
    getCols.push_back("source");
    getCols.push_back("source_fmt");
    getCols.push_back("creation_time");
    getCols.push_back("creator");
    getCols.push_back("quality");
    getCols.push_back("description");
    getCols.push_back("status");

    rdbModel::ResultHandle* res = 0;

    try {
      res = m_rdb->getConnection()->select("Parameters", 
                                           getCols, noCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getParmInfo, key = " << key
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }
    int n = res->getNRows();
    if (!n) {
      delete res;
      return 0;
    }
    
    std::vector<std::string>  fields;
    res->getRow(fields, 0);
    delete res;

    // translate parm class key to parm class name
    std::string parmClassName =
      DbUtil::getColumnValue(m_rdb, "Parameter_class", "name", 
                             "Parameter_class_key", fields[0]);
    return new ParmInfo(keyStr, parmClassName, fields[0], fields[1],
                       fields[2],fields[3],fields[4], fields[5],
                       fields[6],fields[7], fields[8],fields[9]);
  }

  bool MootQuery::getParmsFromMaster(unsigned fmxMasterKey, 
                                     std::vector<ParmOffline>& parms) {
    using facilities::Util;
    parms.clear();
    std::string fmxMasterStr;
    Util::utoa(fmxMasterKey, fmxMasterStr);

    // search in FSW_inputs for  key with FSW_id = fmxMasterKey, status
    // = "added"; return description field. Join with FSW_class to make
    // sure it is a master
    std::string where(" WHERE FSW_id= '");
    where += fmxMasterStr + std::string("' AND status = 'added' ");
    where += ("AND name='LATC_master' AND FSW_class_key=class_fk");

    std::string masterDescrip = 
      DbUtil::getColumnWhere(m_rdb, "FSW_class,FSW_inputs", 
                             "FSW_inputs.description", where, false);

    if (!masterDescrip.size())       return false; // not found

    // masterDescrip should be of form NN_MM_PP_  where NN, MM, etc. are 
    // pos. integers which are themselves FSW_input keys.
    std::vector<std::string> inputKeys;  // for moot keys in FSW_inputs table
    Util::stringTokenize(masterDescrip, std::string("_"), inputKeys);

    for (unsigned ix = 0; ix < inputKeys.size(); ix++) {
      if (!getParmsFromInput(inputKeys[ix], parms)) return false;
    }
    return true;
  }
    
  bool MootQuery::getParmsFromInput(const std::string& fswInputKey, 
                                    std::vector<ParmOffline>& parms) {

    using facilities::Util;

    //  Description field, which will be of the form
    //        pMMpNN..  (i.e., p is delimiter rather than _)
    std::string descrip = 
      DbUtil::getColumnValue(m_rdb, "FSW_inputs", "description", 
                             "FSW_input_key", fswInputKey);

    if (!descrip.size()) return false;
    std::vector<std::string> parmKeys;

    Util::stringTokenize(descrip, std::string("p"), parmKeys);
    unsigned newSize = parms.size() + parmKeys.size();
    parms.reserve(newSize);

    // Make one long WHERE to fetch all the information we need at once
    std::string where("WHERE parm_key = '");
    where += parmKeys[0] + std::string("' ");

    for (unsigned ix = 1; ix < parmKeys.size(); ix++) {
      where += std::string(" OR parm_key = '") + 
        parmKeys[ix] + std::string("' ");
    }
    std::vector<std::string> getCols;
    getCols.reserve(5);
    getCols.push_back("parm_key");
    getCols.push_back("class_fk");
    getCols.push_back("source");
    getCols.push_back("source_fmt");
    getCols.push_back("status");

    rdbModel::ResultHandle* res = 0;
    try {
      res = m_rdb->getConnection()->select("Parameters", getCols, getCols,
                                           where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getParmsFromInput" << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }

    std::string errmsg("MootQuery::getParmsFromInput; wrong # of parameters for fsw key ");
    errmsg += fswInputKey;
    if ((unsigned) DbUtil::checkResults(res, errmsg, parmKeys.size()) 
        != parmKeys.size())
    {
      delete res;
      return false;
    }
    
    for (unsigned i = 0; i < parmKeys.size(); i++) { // stash results
      std::vector<std::string> fields;

      res->getRow(fields, i);
      unsigned pclassNum = Util::stringToUnsigned(fields[1]);

      // Do we need to update cache of parm class names?
      if (pclassNum >= m_parmClasses.size()) updateParmClassCache();
      if (pclassNum >= m_parmClasses.size()) { // TILT!
        // unknown parm class
        delete res;
        return false;  
      }
      // p(key, pclassname, pclassnum, src, srcFmt, status);
      ParmOffline p(fields[0], m_parmClasses[pclassNum], fields[1], 
                    fields[2], fields[3], fields[4]);

      parms.push_back(p);
    }
    delete res;
    return true;
  }

  unsigned MootQuery::getPrecincts(std::vector<std::string>& names) {
    std::string where("");
    names.clear();
    int n = DbUtil::getAllWhere(m_rdb, "Precincts", "name", where,
                                names);
    return (n >= 0) ? 1 : 0;
  }

  VoteAliasInfo* MootQuery::getVoteAliasInfo(unsigned key) {
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    std::string where(" WHERE vote_aliases_key='");
    where += keyStr + std::string("'");

    std::vector<std::string> getCols;
    std::vector<std::string> noCols;
    getCols.reserve(5);
    getCols.clear(); noCols.clear();

    getCols.push_back("name");
    getCols.push_back("precinct_fk");
    getCols.push_back("vote_fk");
    getCols.push_back("last_modified_time");
    getCols.push_back("last_modifier");

    rdbModel::ResultHandle* res = 0;

    try {
      res = m_rdb->getConnection()->select("Vote_aliases", 
                                           getCols, noCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getVoteAliasInfo, key = " << key
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }
    int n = res->getNRows();
    if (!n) {
      delete res;
      return 0;
    }
    
    std::vector<std::string>  fields;
    res->getRow(fields, 0);
    delete res;

    // translate precinct key to precinct name
    std::string precinctName =
      DbUtil::getColumnValue(m_rdb, "Precincts", "name", 
                             "precinct_key", fields[1]);
    return new VoteAliasInfo(keyStr, fields[0], precinctName, fields[1],
                            fields[2],fields[3],fields[4]);
  }

  VoteInfo* MootQuery::getVoteInfo(unsigned key) {
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    std::string where(" WHERE vote_key='");
    where += keyStr + std::string("'");

    std::vector<std::string> getCols;
    std::vector<std::string> noCols;
    getCols.reserve(7);
    getCols.clear(); noCols.clear();

    getCols.push_back("precinct_fk");
    getCols.push_back("instrument");
    getCols.push_back("source");
    getCols.push_back("creation_time");
    getCols.push_back("creator");
    getCols.push_back("description");
    getCols.push_back("status");

    rdbModel::ResultHandle* res = 0;

    try {
      res = m_rdb->getConnection()->select("Votes", getCols, noCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::getVoteInfo, key = " << key
                << " SQL error: "  << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }
    int n = res->getNRows();
    if (!n) {
      delete res;
      return 0;
    }
    
    std::vector<std::string>  fields;
    res->getRow(fields, 0);
    delete res;

    // translate precinct key to precinct name
    std::string precinctName =
      DbUtil::getColumnValue(m_rdb, "Precincts", "name", 
                             "precinct_key", fields[0]);
    return new VoteInfo(keyStr, precinctName, fields[0], fields[1],
                       fields[2],fields[3],fields[4], fields[5],
                        fields[6]);
  }
  bool MootQuery::getVoteParameters(unsigned voteKey, 
                                    std::vector<unsigned>& parmKeys,
                                    bool force) {

    // Get string version of vote key
    std::string voteKeyStr;
    facilities::Util::utoa(voteKey, voteKeyStr);

    parmKeys.clear();
    return voteIsUpToDate(voteKeyStr, &parmKeys, 0, force);
  }

  bool MootQuery::getVotesForPrecinct(const std::string& precinct, 
                                      std::vector<VoteInfo>& voteInfo) {
    // select vote_key from Precincts,Votes where precinct.name=precinct
    std::string where(" WHERE Precincts.name='");
    where += precinct + std::string("'");
    std::vector<unsigned> voteKeys;
    unsigned nKeys;
    try {
      nKeys = DbUtil::getKeys(voteKeys, m_rdb, "Precincts,Votes",
                              "vote_key", where, 0, true);
    }
    catch (std::exception ex) {
      return false;
    }
    voteInfo.reserve(voteInfo.size() + nKeys);
    for (unsigned ix = 0; ix < nKeys; ix++) {
      VoteInfo* v = getVoteInfo(voteKeys[ix]);
      voteInfo.push_back(*v);
      delete v;
    }
    return true;
  }

  unsigned MootQuery::listAncAliasKeys(std::vector<unsigned>& keys,
                                       const std::string& aClass,
                                       const std::string& aliasName) {
    std::string where("");
    static std::string star("*");
    if (aClass != star) {
      where += 
        std::string(" WHERE aclass_fk=(SELECT Ancillary_class_key from Ancillary_class where name='") 
        + aClass + std::string("') ");
    }
    if (aliasName != star) {
      if (where.size() > 0) where += std::string(" AND ");
      else where += std::string(" WHERE ");
      where += std::string(" name='") + aliasName + std::string("'");
    }
    try {
      DbUtil::getKeys(keys, m_rdb, "Ancillary_aliases", 
                      "ancillary_aliases_key", where, 0,
                      false);
    }
    catch (std::exception ex) {
      std::cerr << "Unable to fetch Ancillary alias keys for class=" << aClass
                << ", alias name=" << aliasName << std::endl;
      throw ex;
    }
    return 1;  

  }

  unsigned MootQuery::listAncKeys(std::vector<unsigned>& keys,
                                  const std::string& aClass,
                                  const std::string status,
                                  const std::string quality,
                                  const std::string& instr) {
    std::string where("");
    static std::string star("*");
    if (aClass != star) {
      where += 
        std::string(" WHERE class_fk=(SELECT Ancillary_class_key from Ancillary_class where name='") 
        + aClass + std::string("') ");
    }
    if (status != star) {
      if (where.size() > 0) where += std::string(" AND ");
      else where += std::string(" WHERE ");
      where += std::string(" status='") + status + std::string("'");
    }
    if (quality != star) {
      if (where.size() > 0) where += std::string(" AND ");
      else where += std::string(" WHERE ");
      where += std::string(" quality='") + quality + std::string("'");
    }
    if (instr != star) {
      if (where.size() > 0) where += std::string(" AND ");
      else where += std::string(" WHERE ");
      where += std::string(" instrument='") + instr + std::string("'");
    }
    try {
      DbUtil::getKeys(keys, m_rdb, "Ancillary", "Ancillary_key", where, 0,
                      true);
    }
    catch (std::exception ex) {
      std::cerr << "Unable to fetch Ancillary keys for class=" << aClass
                << ", status=" << status << ", quality=" << quality
                << ", instr=" << instr << std::endl;
      throw ex;
    }
    return 1;  
  }



  unsigned MootQuery::listConfigKeys(std::vector<unsigned>& keys, 
                                     const std::string& status,
                                     const std::string& instr) {
    std::string where("");
    if (status.size()) {
      where = std::string(" WHERE status='" ) + status +
        std::string("'");
    }
    if (instr.size()) {
      if (where.size()) where += std::string(" AND ");
      else where = std::string(" WHERE ");
      where  += std::string("instrument='" ) + instr +
        std::string("'");
    }
    // ask for keys in ascending order
    return DbUtil::getKeys(keys, m_rdb, "Configs", "config_key", where, 
                           0, true);
  }

  unsigned MootQuery::listVoteAliasKeys(std::vector<unsigned>& keys,
                                        const std::string& precinct,
                                        const std::string& aliasName) {
    std::string where("");
    static std::string star("*");
    if (precinct != star) {
      where += 
        std::string(" WHERE precinct_fk=(SELECT precinct_key from Precincts where name='") 
        + precinct + std::string("') ");
    }
    if (aliasName != star) {
      if (where.size() > 0) where += std::string(" AND ");
      else where += std::string(" WHERE ");
      where += std::string(" name='") + aliasName + std::string("'");
    }
    try {
      DbUtil::getKeys(keys, m_rdb, "Vote_aliases", "vote_aliases_key", 
                      where, 0, true);
    }
    catch (std::exception ex) {
      std::cerr << "Unable to fetch Vote_aliases keys for precinct=" 
                << precinct << ", alias name=" << aliasName << std::endl;
      throw ex;
    }
    return 1;  
  }
          


  unsigned MootQuery::listVoteKeys(std::vector<unsigned>& keys,
                                   const std::string& precinct,
                                   const std::string& status,
                                   const std::string& instr) {
    std::string where("");
    static std::string star("*");
    if (precinct != star) {
      where += 
        std::string(" WHERE precinct_fk=(SELECT precinct_key from Precincts where name='") 
        + precinct + std::string("') ");
    }
    if (status != star) {
      if (where.size() > 0) where += std::string(" AND ");
      else where += std::string(" WHERE ");
      where += std::string(" status='") + status + std::string("'");
    }
    try {
      DbUtil::getKeys(keys, m_rdb, "Votes", "vote_key", where, 0,
                      true);
    }
    catch (std::exception ex) {
      std::cerr << "Unable to fetch Vote keys for precinct=" << precinct
                << ", status=" << status 
                << ", instr=" << instr << std::endl;
      throw ex;
    }
    return 1;  
  }



  unsigned MootQuery::resolveAncAlias(const std::string& alias, 
                                      const std::string& ancClass) {
    //                                      ,int tower) {
    // First translate anc class name to key
    std::string where(" WHERE name='");
    where += ancClass + std::string("'");

    std::string ancClassKey = DbUtil::getColumnWhere(m_rdb,
                                                     "Ancillary_class",
                                                     "Ancillary_class_key",
                                                     where, false);

    if (!ancClassKey.size()) return 0;

    return resolveAncAliasByKey(alias, ancClassKey);
  }

  /// Private.  Don't expect external clients to specify by key
  unsigned MootQuery::resolveAncAliasByKey(const std::string& alias, 
                                           const std::string& ancClassKey) {

    std::string where = std::string(" WHERE name='");
    where += alias + std::string("' and aclass_fk='") + ancClassKey
      + std::string("'");
    std::vector<unsigned> keys;  // can only be one returned

    unsigned nKey = DbUtil::getKeys(keys, m_rdb, "Ancillary_aliases",
                                    "ancillary_fk", where);
    return (nKey > 0) ? keys[0] : 0;
  }

  bool MootQuery::resolveAncAliases(std::vector<std::string>& ancKeys,
                                        unsigned voteKey) {
    std::string voteKeyStr;
    facilities::Util::utoa(voteKey, voteKeyStr);
    return resolveAncAliases(ancKeys, voteKeyStr);
  }
  bool MootQuery::resolveAncAliases(std::vector<std::string>& ancKeys,
                                        const std::string& voteKeyStr) {
    std::string where(" WHERE vote_fk='");
    where += voteKeyStr + std::string("' AND aclass_fk IS NOT NULL");

    rdbModel::StringVector getCols;
    getCols.push_back(std::string("aclass_fk"));
    getCols.push_back(std::string("a_alias"));

    rdbModel::ResultHandle* res = 0;
    try {
      res = m_rdb->getConnection()->select("Vote_PClass_AClass", getCols, 
                                           getCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::resolveAncAliases " << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      return -1;
    }

    unsigned nRows = res->getNRows();

    std::vector<std::string> fields;
    fields.reserve(2);    // will hold aclass_fk followed by a_alias

    // Now look each one up in anc alias table
    for (unsigned i = 0; i < nRows; i++) {
      res->getRow(fields, i);
      where = std::string(" WHERE aclass_fk ='") + /*aclassKeys[i]*/ fields[0]
        + std::string("' and name='") 
        + /*aAliases[i]*/ fields[1] + std::string("'");
      std::string ancKey = 
        DbUtil::getColumnWhere(m_rdb, "Ancillary_aliases", "ancillary_fk",
                               where, false);
      // If none found or its blank, failure
      if (ancKey.size() == 0) {
        std::cerr 
          << "MootQuery::resolveAncAliases: Could not resolve all for vote "
          << voteKeyStr << std::endl;
        ancKeys.clear();
        return false;
      }
      ancKeys.push_back(ancKey);
    }
    delete res;
    return true;
  }


  unsigned MootQuery::resolveVoteAlias(const std::string& alias,
                                       const std::string& precinct) {
    // First translate precinct name to key
    std::string where(" WHERE name='");
    where += precinct + std::string("'");

    std::string precinctKey = DbUtil::getColumnWhere(m_rdb,
                                                     "Precincts"
                                                     "precinct_key",
                                                     where, false);

    if (!precinctKey.size()) return 0;

    where = std::string(" WHERE alias='");
    where += alias + std::string("' and precinct_fk='") + precinctKey
      + std::string("'");
    std::vector<unsigned> keys;  // can only be one returned

    unsigned nKey = DbUtil::getKeys(keys, m_rdb, "Vote_aliases",
                                    "vote_aliases_key", where);
    return (nKey > 0) ? keys[0] : 0;
  }


  ////  begin resolveVoteAliases
  bool MootQuery::resolveVoteAliases(std::vector<std::string>& voteKeys,
                                     const std::string& ctnKeyStr) {
    std::string where(" WHERE ctn_fk='");
    where += ctnKeyStr + std::string("' AND precinct_fk IS NOT NULL");

    rdbModel::StringVector getCols;
    getCols.push_back(std::string("precinct_fk"));
    getCols.push_back(std::string("v_alias"));

    rdbModel::ResultHandle* res = 0;
    try {
      res = m_rdb->getConnection()->select("Container_Precinct", getCols, 
                                           getCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::resolveVoteAliases " << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      return -1;
    }

    unsigned nRows = res->getNRows();

    std::vector<std::string> fields;
    fields.reserve(2);    // will hold precinct_fk followed by v_alias
    // Now look each one up in vote alias table
    for (int i = 0; i < nRows; i++) {
      res->getRow(fields, i);

      where = std::string(" WHERE precinct_fk ='") + 
        /* precinctKeys[i] */ fields[0] +
        std::string("' and name='") + /*vAliases[i]*/ fields[1] + 
        std::string("'");
      std::string voteKey = 
        DbUtil::getColumnWhere(m_rdb, "Vote_aliases", "vote_fk",
                               where, false);
      if (voteKey.size() == 0) { // no resolution
        // output complaint?
        std::cerr 
          << "MootQuery::resolveVoteAliases: Could not resolve all for vote "
          << ctnKeyStr << std::endl;
        voteKeys.clear();
        delete res;
        return false;
      }
      voteKeys.push_back(voteKey);
    }
    delete res;
    return true;
  }
  ////  end resolveVoteAliases


  void MootQuery::updateParmClassCache() {
    using facilities::Util;

    std::vector<std::string> getCols;
    getCols.push_back("Parameter_class_key");
    getCols.push_back("name");

    rdbModel::ResultHandle* res = 0;
    try {
      res = m_rdb->getConnection()->select("Parameter_class", getCols,
                                           getCols);
    }
    catch (std::exception ex) {
      std::cerr << "MootQuery::updateParmClassCache SQL error: " 
                << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
    }
    
    m_parmClasses.clear();

    // find value of last key
    // resize m_parmClasses accordingly.  Then store
    // name values in vector spot corresponding to key number
    unsigned nrows = res->getNRows();
    std::vector<std::string> fields;
    res->getRow(fields, nrows - 1);
    unsigned lastKey = Util::stringToUnsigned(fields[0]);
    m_parmClasses.resize(lastKey + 1);
    
    for (unsigned ix = 0; ix < nrows; ix++) {
      res->getRow(fields, ix);
      unsigned key = Util::stringToUnsigned(fields[0]);
      m_parmClasses[key] = fields[1];
    }

    if (res) delete res;
  }

  bool MootQuery::voteExists(unsigned voteKey, bool goodStatus) {
    std::string voteKeyStr;
    facilities::Util::utoa(voteKey, voteKeyStr);
    return voteExists(voteKeyStr, goodStatus);
  }

  bool MootQuery::voteExists(const std::string& voteKeyStr, bool goodStatus) {

    std::string voteStatus;
    try {
      voteStatus = 
        DbUtil::getColumnValue(m_rdb, "Votes", "status", "vote_key", 
                               voteKeyStr); 
    }
    catch (DbUtilNoDataException ex) {
      std::cerr << "MootQuery::voteExists:  " << voteKeyStr 
                << " is not the key of a valid registered vote" << std::endl;
      return false;
    }
    if (goodStatus) {  // status must be "CREATED"
      if (voteStatus != std::string("CREATED") ) {
        std::cerr << "MootQuery::voteIsContainer:  " << voteKeyStr 
                  << " is not the key of a valid registered vote" << std::endl;
        return false;
      }
    }
    return true;
  }

  bool MootQuery::voteIsUpToDate(unsigned voteKey, std::vector<unsigned>* pk,
                                 bool* isCtn, bool force) {

    std::string voteKeyStr;
    facilities::Util::utoa(voteKey, voteKeyStr);

    return voteIsUpToDate(voteKeyStr, pk, isCtn, force);
  }

  /**
     A vote is up to date if
        - for each parameter class controlled by the vote, there is
          a registered parameter file instance (row in Parameters table)
          of suitable quality which references this vote with its 
          vote_fk field
        - if the vote file specifies that a parameter class depends on
          one or more ancillary classes, there must be an instance of
          the class 
           * which references the vote with its vote_fk field, as above
           * whose key appears in  Parameters_to_Ancillary exactly once
             for each ancillary class it should depend on, according
             to the vote file.  The corresponding Ancillary_fk entries
             in those rows must refer to the same ancillary files one
             would get by resolving the ancillary aliases mentioned
             in the vote file

             If @a parmKeys is non-zero, fill the vector with the
             parameter files (one per class) associated with the vote.
   */
  bool MootQuery::voteIsUpToDate(const std::string& voteKeyStr,
                                 std::vector<unsigned>* pk, bool* isCtn,
                                 bool force) {

    // First be sure vote is registered and valid
    bool isContainer;
    bool ret = true;

    if (!voteExists(voteKeyStr, &isContainer) ) {
      std::cerr << "MootQuery::voteIsUpToDate:  " << voteKeyStr 
                << " is not the key of a valid registered vote" << std::endl;
      return false;
    }

    if (pk) pk->clear();

    // First find out whether we're dealing with a container or not
    std::vector<unsigned> precinctKeys;
    std::string where(" WHERE ctn_fk ='");
    where += voteKeyStr + std::string("'");
    unsigned nContained = 
      DbUtil::getKeys(precinctKeys, m_rdb, "Container_Precinct",
                      "precinct_fk", where);
    if (isCtn) *isCtn = (nContained > 0);   // inform caller 

    if (nContained) {
      std::vector<std::string> innerKeys;
      std::vector<unsigned> parmKeys;

      if (!resolveVoteAliases(innerKeys, voteKeyStr)) return false;
      for (unsigned iVote = 0; iVote < nContained; iVote++) {
        bool goodPrecinct = 
          voteIsUpToDate(innerKeys[iVote], &parmKeys, 0, force);
        //        if (voteIsUpToDate(innerKeys[iVote], &parmKeys, 0, force)) {

        if (!goodPrecinct) {
          if (!parmKeys.size()) {   // Done. This vote was never even built
            std::cerr << "Container vote " << voteKeyStr 
                      << "not up to date; failed for unbuilt vote "
                      << innerKeys[iVote] << std::endl;
            if (pk) pk->clear();
            return false;
          }
          std::cerr << "Container vote " << voteKeyStr 
                    << "not up to date; failed for resolved vote "
                    << innerKeys[iVote] << std::endl;
          ret = false;
        }
        if (goodPrecinct || force) { 
          if (pk) pk->insert(pk->end(), parmKeys.begin(), parmKeys.end());
        }
        else {  // bad precinct, no force.  We're done
          if (pk) pk->clear();
          return false;
        }
      }
      return ret;
    }

    std::vector<std::string> aliasedAncKeys;

    static std::string goodParm("' AND status='CREATED' AND quality='PROD'");

    if (!resolveAncAliases(aliasedAncKeys, voteKeyStr) ) {
      return false;
    }

    // Get parameter classes associated with this vote.
    std::vector<std::string> pclassKeyStr;


    where = std::string(" WHERE vote_fk ='");
    where += voteKeyStr + std::string("' AND aclass_fk IS NULL");
    int nPclass = DbUtil::getAllWhere(m_rdb, "Vote_PClass_AClass",
                                      "pclass_fk", where, pclassKeyStr);

    if (!nPclass) {
      std::cerr << "Vote file "  << voteKeyStr << 
        " has no associated parameter classes" << std::endl;
      return false;
    }
    // Find Parameter file *instances* associated with this vote
    // Must be at least one per class.
    // What about instrument?  Should we be checking that too?
    // If so caller will have to pass it in.
    // Currently don't have one associated with a MootQuery object.
    for (int iPclass=0; iPclass < nPclass; iPclass++) {
      where = std::string(" WHERE vote_fk ='") + voteKeyStr +
        std::string("' and class_fk='") + pclassKeyStr[iPclass] + goodParm;

      std::vector<std::string> parmKeys;
      int nParm = 
        DbUtil::getAllWhere(m_rdb, "Parameters", "parm_key", where, parmKeys);
      if (nParm < 0) throw std::runtime_error("VoteIsUpToDate db error");
      else if (nParm == 0) return false;

      std::vector<std::string> ancClasses;
      where = std::string(" WHERE vote_fk ='") + voteKeyStr +
        std::string("' and pclass_fk='") + pclassKeyStr[iPclass] + 
        std::string("' AND aclass_fk IS NOT NULL");
      
      //    From Vote_PClass_AClass retrieve anc classes the pclass depends
      //    on.  If none, we're done (successfully) with this parm class.
      int nAncClass = DbUtil::getAllWhere(m_rdb, "Vote_PClass_AClass",
                                     "aclass_fk", where, ancClasses);
      if (nAncClass < 0) { 
        throw std::runtime_error("VoteIsUpToDate db error");
      }
      else if (nAncClass == 0) {
        if (pk) {
          unsigned lastKey = 
            facilities::Util::stringToUnsigned(parmKeys.back());
          pk->push_back(lastKey);
        }
        continue; // no dependencies for this parm
      }
      bool parmClassOk = false;

      // Otherwise, consider all "good" instances of parm class associated
      // with this vote, fetched above.  More likely to succeed with 
      // most recent, so start there
      unsigned iP = parmKeys.size();
      while (iP--) {               // look at collection of parameter instances
        std::string pkey = parmKeys[iP];
        bool goodParmInstance = true;
        //  fetch anc keys from rows in Parameters_to_Ancillary with
        //  parm key = the one we're examining.
        std::vector<std::string> inputAncKeys;
        where = std::string(" WHERE Parameter_fk='") + pkey + std::string("'");
        int nAnc = DbUtil::getAllWhere(m_rdb, "Parameters_to_Ancillary",
                                       "Ancillary_fk", where, inputAncKeys);
        // Param should be associated with nAncClass anc files *and* they
        // should be of the correct classes as specified in the vote file
        // (equivalently, in table Vote_PClass_AClass)
        // **TODO** add this check.    So far, just checking that the
        // count is correct.  nAncClass is the number of ancillaries the
        // vote file says the parameter should depend on; nAnc is the 
        // number this particular parameter file does depend on.
        if (nAnc != nAncClass) {
          goodParmInstance = false;
          continue; 
        }
        for (int iAnc = 0; iAnc < nAnc; iAnc++) {
          if (find(aliasedAncKeys.begin(), aliasedAncKeys.end(), 
                   inputAncKeys[iAnc]) 
              == aliasedAncKeys.end() ) {
            goodParmInstance = false;
            break;   // try another 
          }
        }
        if (goodParmInstance) {
          parmClassOk = true;
          if (pk) {   // save for caller
            unsigned lastKey = 
              facilities::Util::stringToUnsigned(pkey);
            pk->push_back(lastKey);
          }
          break;
        }
      //      If all anc keys are among set of resolved aliases, we're
      //      done (successfully) with this pclass
      }     // end 'look' at a collection of instances for a param. class
      if (!parmClassOk) {
        if (force) {  // keep going, but will ultimately return bad status
          ret = false;
        }
        else {
          if (pk) pk->clear();
          return false;
        }
      }
    }  // got through all parm classes one way or another
    return ret;
  }

}
