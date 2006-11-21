//  $Header$

#include <string>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "mootCore/MootQuery.h"
#include "mootCore/MoodConnection.h"
#include "mootCore/DbUtil.h"
#include "rdbModel/Rdb.h"
#include "rdbModel/Db/Connection.h"
#include "rdbModel/Db/ResultHandle.h"
#include "facilities/Util.h"
#include "mootCore/DbUtilException.h"

namespace MOOT {
  MootQuery::MootQuery(MoodConnection* mood) : m_rdb(0), m_mood(mood), 
                                               m_ownConnection(false),
                                               m_archive("")  {
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

  bool MootQuery::getConfigFswInfoByAlg(const std::string& algName,
                                        unsigned algStep,
                                        std::vector<FswInfo>& info,
                                        bool clear) {
    std::string stepStr;
    facilities::Util::utoa(algStep, stepStr);
    return getConfigFswInfoByAlg(algName, stepStr, info, clear);
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
 // FswInfo(fmxPath,src,description,fmxDb,fmxKey,mootClassKey);
      FswInfo f(selFields[0], selFields[1], selFields[2], selFields[5],
                id, classFk);
      info.push_back(f);
      delete res;
    }
    return true;
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

  std::string MootQuery::classStr(const std::string& table, unsigned key) {
    std::string keyStr;
    facilities::Util::utoa(key, keyStr);
    std::string keyCol = table + "_key";
    try {
      return DbUtil::getColumnValue(m_rdb, table, "name", keyCol, keyStr);
    }
    catch (std::exception ex) {
      std::cerr << ex.what() << std::endl;
      std::cerr.flush();
      throw ex;
    }
  }

  int MootQuery::classKey(const std::string& table, const std::string& name) {
    std::string getCol = table + "_key";
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
    getCols.reserve(9);
    getCols.push_back("config_key");
    getCols.push_back("name");
    getCols.push_back("algorithm");
    getCols.push_back("alg_step");
    getCols.push_back("description");
    getCols.push_back("status");
    getCols.push_back("active_state");
    getCols.push_back("mode");
    getCols.push_back("creation_request_time");

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
                       selFields[8]);
      info.push_back(info1);
    }
    return n;
  }
}
