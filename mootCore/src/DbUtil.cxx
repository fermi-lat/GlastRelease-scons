// $Header$

#include <cstdio>
// #include <cstdlib>
#include <iostream>
#include "facilities/Util.h"
#include "rdbModel/Rdb.h"
#include "rdbModel/Db/MysqlConnection.h"
#include "rdbModel/Db/ResultHandle.h"
#include "rdbModel/Tables/Table.h"
#include "rdbModel/Tables/Column.h"
#include "rdbModel/Tables/Assertion.h"
#include "rdbModel/RdbException.h"
#include "mootCore/DbUtil.h"
#include "mootCore/DbUtilException.h"

namespace {
  bool connected(rdbModel::Rdb* rdb) {
    if (!rdb) {
      std::cerr << "DbUtil:  need Rdb object" 
                << std::endl;
      std::cerr.flush();
      return false;
    }
    if (!rdb->getConnection()) {
      std::cerr << "DbUtil: need connection to db"
                << std::endl;
      std::cerr.flush();
      return false;
    }
    return true;
  }
}

namespace MOOT {
  /** 
      Check that SELECT completed ok. If so return # rows. Else return -1
      If 'expected' > 0, check #rows actually returned = expected.
      If 'expected' < 0, check # rows returned >= (-expected)
      If 'expected = 0, don't check.
  */
  int DbUtil::checkResults(rdbModel::ResultHandle* res, 
                           const std::string& errMsgPrefix, int expected) {
    if (!res) {
      return -1;
      //      std::cerr << errMsgPrefix << " bad query" << std::endl;
      //      std::cerr.flush();
      //      return 0;
    }
    int got = res->getNRows();
    if (expected > 0 ) {
      if (expected != got) {
        std::cerr << errMsgPrefix << " expected " << expected 
                  << " rows; got " << got << std::endl;
        std::cerr.flush();
      }
    }
    else if (expected < 0) {
      expected = -expected;
      if (got < expected) {
        std::cerr << errMsgPrefix << " expected at least" << expected 
                  << " rows; got " << got << std::endl;
        std::cerr.flush();
      }
    }
    return got;
  }

  std::string DbUtil::getColumnValue(rdbModel::Rdb* rdb,
                                     const std::string& table, 
                                     const std::string& col,
                                     const std::string& keyCol,
                                     const std::string& key,
                                     bool forUpdate) {

    if (!connected(rdb)) {
      throw rdbModel::RdbException("DbUtil::getColumnValue: bad connection");
    }
      //return std::string("");
    std::string where(" WHERE ");
    where += keyCol + std::string("= '") + key + std::string("'");
    return getColumnWhere(rdb, table, col, where, true, forUpdate);
  }

  std::string DbUtil::getColumnWhere(rdbModel::Rdb* rdb,
                                     const std::string& table, 
                                     const std::string& col,
                                     const std::string& where, 
                                     bool onlyOne,
                                     bool forUpdate) {

    if (!connected(rdb))  {
      throw rdbModel::RdbException("DbUtil::getColumnValue: bad connection");
    }

      // return std::string("");

    std::string myWhere(where);
    rdbModel::StringVector getCols;
    getCols.push_back(col);
    rdbModel::StringVector noCols;
    bool useNoCols = false;
    if (forUpdate) {
      myWhere += " FOR UPDATE ";
      useNoCols = true;
    }
    if (myWhere.find("ORDER BY") < myWhere.size()) useNoCols = true;
    rdbModel::ResultHandle* res = 0;
    try {
      if (useNoCols) {// prevent MysqlConnection::select from adding ORDER BY
        res = rdb->getConnection()->select(table, getCols, noCols, myWhere);
      } else {
        res = rdb->getConnection()->select(table, getCols, getCols, myWhere);
      }
    }
    catch (std::exception ex) {
      std::cerr << "DbUtil::getColumnWhere " << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      throw ex;
      //      return std::string("");
    }
    int nRows = checkResults(res, "MootBuild,getColumnValue ", 0);
    if (nRows == 0) {
      if (!onlyOne) return std::string("");
      else throw DbUtilNoDataException("DbUtil::getcolumnWhere: no data ");
    }
    else if (nRows < 0) 
    {
      throw DbUtilException("DbUtil::getColumnWhere:  query failed ");
    }
    if ((nRows > 1) && onlyOne ) {
      std::cerr << "DbUtil::getColumnWhere: too many rows satisfy "
                << myWhere << std::endl;
      std::cerr.flush();
      delete res;
      //      return std::string("");
      throw DbUtilNotUniqueException("DbUtil::getColumnWhere: too many found");
    }
    std::vector<std::string>selFields;

    res->getRow(selFields);
    std::string val = selFields[0];
    delete res;
    return val;
  }

  int DbUtil::getAllWhere(rdbModel::Rdb* rdb,
                          const std::string& table,
                          const std::string& col,
                          const std::string& where,
                          std::vector<std::string>& vals,
                          bool atLeastOne) {
    if (!connected(rdb)) return -1;

    rdbModel::StringVector getCols;
    getCols.push_back(col);

    rdbModel::ResultHandle* res = 0;
    try {
      res = rdb->getConnection()->select(table, getCols, 
                                         getCols, where);
    }
    catch (std::exception ex) {
      std::cerr << "MootBuild::getAllWhere " << ex.what() << std::endl;
      std::cerr.flush();
      if (res) delete res;
      return -1;
    }

    unsigned nRows = res->getNRows();
    int nFetched = 0;

    for (unsigned ix = 0; ix < nRows; ix++) {
      if (res->getRow(vals, ix, false)) {
        nFetched++;
      }
    }
    delete res;

    if (( res <= 0) && atLeastOne)
      throw DbUtilException("DbUtil::getColumnValue: no valid match");

    return nFetched;
  }

  unsigned DbUtil::getKeys(std::vector<unsigned>& keys,
                           rdbModel::Rdb* rdb, const std::string& table,
                           const std::string& keyCol, 
                           const std::string& where, unsigned limit,
                           bool ascend) {
    std::vector<std::string> getCols;
    getCols.push_back(keyCol);
    rdbModel::SELECTOPTIONS ordering = (ascend) ? rdbModel::SELECTnone 
      : rdbModel::SELECTdesc;
    rdbModel::ResultHandle* res = 
      rdb->getConnection()->select(table, getCols, getCols,
                                   ordering, where, limit);

    unsigned nRows = res->getNRows();
    std::vector<std::string> fields;
    for (unsigned iRow = 0; iRow < nRows; iRow++) {
      res->getRow(fields, iRow);
      keys.push_back(facilities::Util::stringToUnsigned(fields[0]));
    }
    delete res;
    return nRows;
  }

  int DbUtil::setColumnValue(rdbModel::Rdb* rdb,
                            const std::string& table,
                            const std::string& colName,
                            const std::string& value,
                             const std::string& where) {
    using rdbModel::FieldVal;
    using rdbModel::Row;

    std::vector<FieldVal> fields;
    fields.push_back(FieldVal(colName, value, ""));
    Row row(fields);
    return rdb->updateRows(table, row, where);
  }

  int DbUtil::setColumnValue(rdbModel::Rdb* rdb,
                             const std::string& table,
                             const std::string& colName,
                             const std::string& value,
                             const std::string& keyCol,
                             const std::string& key) {
    std::string where(" WHERE ");
    where += keyCol + std::string("= '") + key + std::string("'");
    return setColumnValue(rdb, table, colName, value, where);
  }

  bool DbUtil::startTransaction(rdbModel::Rdb* rdb) {
    try {
      rdb->getConnection()->dbRequest("start transaction");
    }
    catch (...) {
      std::cerr << "DbUtil::startTransaction: request failed" << std::endl;
      return false;
    }
    return true;
  }

  bool DbUtil::commit(rdbModel::Rdb* rdb) {
    try {
      rdb->getConnection()->dbRequest("commit");
    }
    catch (...) {
      std::cerr << "DbUtil::commit: request failed" << std::endl;
      return false;
    }
    return true;
  }

}
