// $Header$
#ifdef  WIN32
#include <windows.h>
#endif

#include "rdbModel/Db/MysqlConnection.h"
#include "rdbModel/Rdb.h"
#include "rdbModel/Tables/Table.h"
#include "rdbModel/Tables/Assertion.h"
#include "rdbModel/Tables/Column.h"
#include "rdbModel/Tables/Datatype.h"
#include "rdbModel/Db/MysqlResults.h"
#include "rdbModel/RdbException.h"
#include "facilities/Util.h"

#include "xml/XmlParser.h"
#include "xml/Dom.h"

// #include "rdbModel/Management/Manager.h"

#include "mysql.h"
#include <iostream>
#include "facilities/Util.h"
namespace {

  // Size specification is of form (m) or (m,d)  If no size specification 
  // return 0; else return value of m.
  int extractSize(const std::string& sqlString) {
    unsigned leftLoc = sqlString.find("(");
    if (leftLoc == std::string::npos) return 0;
    leftLoc++;           // now is at start of m
    unsigned rightLoc = sqlString.find(",");
    if (rightLoc == std::string::npos) {
      rightLoc = sqlString.find(")");
    }
    std::string numString = 
      sqlString.substr(leftLoc, rightLoc - leftLoc);
    return facilities::Util::stringToInt(numString);
  }

  void addArg(bool literal, const std::string arg, std::string& sqlString) {
    if (literal) sqlString += '"';
    sqlString += arg;
    if (literal) sqlString += '"';
    return;
  }

  bool compareEnumList(const std::vector<std::string>& choices, 
                       std::string sqlType) {
    // Number has to be the same.  
    unsigned locComma = sqlType.find(",");
    unsigned nComma = 0;
    while (locComma != std::string::npos) {
      nComma++;
      locComma = sqlType.find(",", locComma+1);
    }
    unsigned nChoice = choices.size();
    if (nChoice != (nComma + 1)) return false;
    for (unsigned iChoice = 0; iChoice < nChoice; iChoice++) {
      unsigned loc = sqlType.find(choices[iChoice]);
      if (loc == std::string::npos) return false;
    }
    return true;
  }
}
    
namespace rdbModel {
  bool   MysqlConnection::m_compileInit = false;

  MysqlConnection::MysqlConnection(std::ostream* out,
                                   std::ostream* errOut) :
    m_mysql(0), m_connected(0), m_out(out), m_err(errOut),
    m_visitorType(VISITORundefined), m_tempRes(0) {
    if (m_out == 0) m_out = &std::cout;
    if (m_err == 0) m_err = &std::cerr;
  }

  bool MysqlConnection::close() {
    if (m_tempRes) {
      mysql_free_result(m_tempRes);
      m_tempRes = 0;
    }
    mysql_close(m_mysql);
    m_mysql = 0;
    m_connected = false;
    return true;
  }

  MysqlConnection::~MysqlConnection() {
    close();
    delete m_mysql;
    return;
  }

  bool MysqlConnection::open(const std::string& host, 
                             const std::string& user,
                             const std::string& password,
                             const std::string& dbName,
                             unsigned int       port) {
    if (dbName.size() == 0) {
      (*m_err) << 
        "rdbModel::MysqlConnection::open : null db name not allowed!" <<
        std::endl;
      return false;
    } 

    m_mysql = new MYSQL;
    mysql_init(m_mysql);

    //    mysql_init(m_mysql);
    MYSQL *connected = mysql_real_connect(m_mysql, host.c_str(), user.c_str(),
                                          password.c_str(), dbName.c_str(),
                                          port, NULL, 0);

    if (connected != 0) {  // Everything is fine.  Put out an info message
      (*m_out) << "Successfully connected to MySQL host " << 
        host << std::endl;
      m_connected = true;
      m_dbName = dbName;
    }
    else {
      (*m_err) <<  "Failed to connect to MySQL host " << host <<
        "with error " << mysql_error(m_mysql) << std::endl;
      m_connected = false;
    }
    return m_connected;
  }

  bool MysqlConnection::open(const std::string& parms) {

    xml::XmlParser parser;
    DOM_Document doc = parser.parse(parms.c_str(), "mysqlConnection");
    if (doc == DOM_Document()) {
      (*m_err) << "parse of connection parameters failed" << std::endl;
      return false;
    }
    DOM_Element  conn = doc.getDocumentElement();
    
    std::string host = xml::Dom::getAttribute(conn, "host");
    std::string user = xml::Dom::getAttribute(conn, "user");
    std::string password = xml::Dom::getAttribute(conn, "password");
    std::string dbname = xml::Dom::getAttribute(conn, "dbname");
    int port = xml::Dom::getIntAttribute(conn, "port");
    if (password.size() == 0 ) { // prompt for password?
      (*m_out) << "interactive login NYI " << std::endl;
      return false;
    }
    return this->open(host, user, password, dbname, port);
  }



  MATCH MysqlConnection::matchSchema(Rdb *rdb) {
    if (!m_connected) return MATCHnoConnection;

    // Check global characteristics; 
    // Could do this via Manager; seems a bit artificial, bypass for now
    m_visitorType = VISITORmatch;
    m_matchReturn = MATCHequivalent;
    unsigned int ret = rdb->accept(this);

    if ((ret == Visitor::VERROR) || (ret == Visitor::VERRORABORT)) {
      return MATCHfail;
    }
    else return m_matchReturn;
  }



    // For each table
    //         compare # of columns
    //         compare datatype description, other attributes of column
    //         compare indices

  



  bool MysqlConnection::insertRow(const std::string& tableName, 
                                  const StringVector& colNames, 
                                  const StringVector& values,
                                  int* auto_value) {
    std::string ins;
    if (auto_value) *auto_value = 0;

    // check that sizes of vectors match
    unsigned  nCol = colNames.size();    
    if (!nCol || (nCol != values.size()  ) ) {
      (*m_err) << " MysqlConnection::insertRow: vector lengths incompatible"
                << std::endl;
      return false;
    }

    // caller should already have checked for validity and should
    // have supplied all necessary columns

    ins += "insert into " + tableName;
    ins += "set " + colNames[0] + "='" + values[0] + "' ";
    for (unsigned iCol = 1; iCol < nCol; iCol++) {
      ins += ",set " + colNames[iCol] + "='" + values[iCol] + "' ";
    }

    int mysqlRet = mysql_query(m_mysql, ins.c_str());

    if (mysqlRet) {
      (*m_err) << "MySQL error during INSERT, code " << mysqlRet << std::endl;
      return false;
    }
    if (auto_value) {
      *auto_value = mysql_insert_id(m_mysql);
    }
    return true;
  }


  unsigned int MysqlConnection::update(const std::string& tableName, 
                                       const StringVector& colNames, 
                                       const StringVector& values,
                                       const Assertion* where) {

    unsigned int nCol = colNames.size();
    if (nCol != values.size()) {
      (*m_err) << "rdbModel::mysqlConnection::update: ";
      (*m_err) << "Incompatible vector arguments " << std::endl;
      return 0;
    }
    std::string sqlString = "UPDATE " + tableName + " SET ";
    sqlString += colNames[0] + "'" + values[0] + "'";
    for (unsigned int iCol = 1; iCol < nCol; iCol++) {
      sqlString += "," + colNames[iCol] + "'" + values[iCol] + "'";
    }
    if (where) {
      sqlString += " WHERE ";
      bool ret = compileAssertion(where, sqlString);
      if (!ret) return 0;
    }
    int mysqlRet = mysql_query(m_mysql, sqlString.c_str());

    if (mysqlRet) {
      (*m_err) << "rdbModel::MysqlConnection::update: ";
      (*m_err) << "MySQL error during UPDATE, code " << mysqlRet << std::endl;
      return 0;
    }
    my_ulonglong nModLong = mysql_affected_rows(m_mysql);
    // Not much chance that we'll change more rows than will fit in just long
    unsigned nMod = nModLong;
    return nMod;


  }

  ResultHandle* MysqlConnection::select(const std::string& tableName,
                                        const StringVector& getCols,
                                        const StringVector& orderCols,
                                        const Assertion* where,
                                        int   rowLimit) {
    std::string sqlString = "SELECT FROM " + tableName +  " ";
    unsigned nGet = getCols.size();
    unsigned nOrder = orderCols.size();

    sqlString += getCols[0];
    for (unsigned iGet = 1; iGet < nGet; iGet++) {
      sqlString += ",";
      sqlString += getCols[iGet];
    }
    if (where != 0) {
      sqlString += " WHERE ";
      bool ret = compileAssertion(where, sqlString);
      if (!ret) return 0;
    }
    if (nOrder > 0 ) {
      sqlString += " ORDER BY " + orderCols[0]; 
      for (unsigned iOrder = 1; iOrder < nOrder; iOrder++) {
        sqlString += ",";
        sqlString += orderCols[iOrder];
      }
    }
    if (rowLimit > 0) {
      std::string limitStr;
      facilities::Util::itoa(rowLimit, limitStr);
      sqlString += " LIMIT " + limitStr;
    }
    int mysqlRet = mysql_query(m_mysql, sqlString.c_str());
    if (mysqlRet) {
      std::string msg = 
        "rdbModel::MysqlConnection::select: mysql_query error, code ";
      std::string codeString;
      facilities::Util::itoa(mysqlRet, codeString);
      msg += codeString;
      throw RdbException(msg, mysqlRet);
      return 0;
    }

    MYSQL_RES *myres = mysql_store_result(m_mysql);
    MysqlResults* results = new MysqlResults(myres);
    return results;
  }

  ResultHandle* MysqlConnection::dbRequest(const std::string& request) {
    int mysqlRet = mysql_query(m_mysql, request.c_str());
    if (mysqlRet) {
      std::string msg = 
        "rdbModel::MysqlConnection::dbRequest: mysql_query error, code ";
      std::string codeString;
      facilities::Util::itoa(mysqlRet, codeString);
      msg += codeString;
      throw RdbException(msg, mysqlRet);
      return 0;
    }

    MYSQL_RES *myres = mysql_store_result(m_mysql);
    if (!myres) {
      // Was it supposed to return data?
      if (mysql_field_count(m_mysql) == 0) { // no data expected
        return 0;
      }
      else {
        std::string msg =
          "rdbModel::MysqlConnection::dbRequest: expected data; none returned";
        throw RdbException(msg);
        return 0;
      }
    }
    return new MysqlResults(myres);
  }

  bool MysqlConnection::compileAssertion(const Assertion* a, 
                                         std::string& sqlString) const {
    if (!m_compileInit) {
      compileInit();
      m_compileInit = true;
    }
    return compileOperator(a->getOperator(), sqlString);
  }

  std::string opSymbols[OPTYPElast];

  void MysqlConnection::compileInit() {
    opSymbols[OPTYPEor] = " OR ";
    opSymbols[OPTYPEand] = " AND ";
    opSymbols[OPTYPEnot] = " NOT ";
    opSymbols[OPTYPEexists] = "EXISTS ";
    opSymbols[OPTYPEisNull] = " IS NULL";
    opSymbols[OPTYPEequal] = "=";
    opSymbols[OPTYPEnotEqual] = "<>";
    opSymbols[OPTYPElessThan] = "<";
    opSymbols[OPTYPEgreaterThan] = ">";
    opSymbols[OPTYPElessOrEqual] = "<=";
    opSymbols[OPTYPEgreaterOrEqual] = ">=";
    return;
  }

  /** Result is appended to caller-supplied string
   Convention is to use "   "  around literal values
   Note no verification is done here; that operator is in fact a comparison
   or isNull operator.  This is called internally only and that check will
   have been done before invoking this routine.
  */
  bool MysqlConnection::compileComparison(Assertion::Operator* op, 
                                          std::string& sqlString) {
    OPTYPE opType = op->getOpType();
    if (opType == OPTYPEisNull) {
      sqlString +="(";
      sqlString += op->getCompareArgs()[0];
      sqlString += opSymbols[opType];
      sqlString += ")";
      return true;
    }
    sqlString += "(";
    addArg(op->getLiteralness()[0], op->getCompareArgs()[0], sqlString);
    sqlString += opSymbols[opType];
    addArg(op->getLiteralness()[1], op->getCompareArgs()[1], sqlString);
    sqlString += ")";

    return true;
  }

  bool MysqlConnection::compileOperator(Assertion::Operator* op, 
                                        std::string &sqlString) {
    if (op->isCompareOp() ) return compileComparison(op, sqlString);
    
    bool ret = true;

    const std::vector<Assertion::Operator*>& children = op->getChildren();
    unsigned nChild = children.size();

    sqlString += "(";
    // For single-child operators NOT,  exists, operator symbol
    // goes 1st, then operand
    if (nChild <= 1) { // operator goes first
      sqlString += opSymbols[op->getOpType()];

      // more special handling for EXISTS
      if (op->getOpType() == OPTYPEexists) {
        sqlString += "(SELECT * FROM " + op->getTableName();
        if (!nChild) {     // done
          sqlString += ")";
          return ret;
        }
        // else EXISTS child is object of a WHERE clause 
        sqlString += " WHERE(";
      }
      ret = compileOperator(children[0], sqlString);
      if (!ret) {
        std::string msg = 
          "rdbModel::MysqlConnection::compileOperator failed for operator "
          + opSymbols[op->getOpType()];
        throw RdbException(msg);
      }
      sqlString += ")";

      // Have an extra closing ")" for EXISTS with WHERE clause
      if (op->getOpType() == OPTYPEexists)       sqlString += ")";

      return ret;
    }
    // Otherwise put operator symbols between adjacent children.
    std::string symbol = opSymbols[op->getOpType()];

    ret = compileOperator(children[0], sqlString);
    if (!ret) {
      std::string msg = 
        "rdbModel::MysqlConnection::compileOperator failed for operator "
        + symbol;
      throw RdbException(msg);
    }
    for (unsigned int iChild = 1; iChild < nChild; iChild++) {
      sqlString += symbol;
      ret = compileOperator(children[iChild], sqlString);
      if (!ret) {
        std::string msg = 
          "rdbModel::MysqlConnection::compileOperator failed for operator "
          + symbol;
        throw RdbException(msg);
      }
    }
    return ret;
  }

  // Satisfy Visitor interface.  For now the only visitor is the
  // one to check whether remote and local db descriptions match or
  // are at least compatible enough to be used.
  Visitor::VisitorState MysqlConnection::visitRdb(Rdb *rdb) {
    
    if (m_dbName != rdb->getDbName()) {
      m_matchReturn = MATCHfail;
      return Visitor::VDONE;
    }
    
    unsigned int nLocal = rdb->getNTable();

    // Null pointer for 2nd argument means "list all tables"
    
    MYSQL_RES* res = mysql_list_tables(m_mysql, 0);
    if (!res) {
      m_matchReturn = MATCHfail;
      return Visitor::VERRORABORT;
    }
    unsigned int nRemote = mysql_num_rows(res);
    mysql_free_result(res);

    if (nRemote < nLocal) {
      m_matchReturn = MATCHfail;
      return Visitor::VDONE;
    }
    else if (nRemote > nLocal) m_matchReturn = MATCHcompatible;

    return Visitor::VCONTINUE;
  }

  Visitor::VisitorState MysqlConnection::visitTable(Table* table) {
    const std::string& tName = table->getName();

    // Result set will have all fields for the table
    if (m_tempRes) {
      mysql_free_result(m_tempRes);
      m_tempRes = 0;
    }
    m_primColName.clear();

    std::string query = "SHOW COLUMNS FROM " + tName;

    int ret = mysql_query(m_mysql, query.c_str());
    if (ret) {
      m_matchReturn = MATCHfail;
      return Visitor::VERRORABORT;
    }
      
    m_tempRes = mysql_store_result(m_mysql);
    if (!m_tempRes) {
      m_matchReturn = MATCHfail;
      return Visitor::VERRORABORT;
    }
    // Result set is a table with fields "Field"(the name) "Type" "Null"(yes
    // or no) "Key" "Default", "Extra"  
    // Make it easier for accept(Column* ) to find relevant information
    unsigned int nRow = mysql_num_rows(m_tempRes);
    m_colIx.clear();
    for (unsigned iRow = 0; iRow < nRow; iRow++) {
      MYSQL_ROW colDescrip = mysql_fetch_row(m_tempRes);
      std::string name = std::string(colDescrip[0]);
      m_colIx[name] = iRow;
    }
    return Visitor::VCONTINUE;

  }

  Visitor::VisitorState MysqlConnection::visitColumn(Column* col) {
    std::string myName = col->getName();
    if (m_colIx.find(myName) == m_colIx.end()) {
      m_matchReturn = MATCHfail;
      return Visitor::VERRORABORT;
    }
    unsigned int ix = m_colIx[myName];
    mysql_data_seek(m_tempRes, ix);
    MYSQL_ROW colDescrip = mysql_fetch_row(m_tempRes);

    // Type
    std::string sqlDtype = std::string(colDescrip[1]);
    Datatype* dtype = col->getDatatype();
    if (!checkDType(dtype, sqlDtype)) {
      m_matchReturn = MATCHfail;
      return Visitor::VERRORABORT;
    }
    
    // Null
    bool nullable = (std::string(colDescrip[2]) == std::string("YES"));
    if (nullable != col->nullAllowed()) {
      m_matchReturn = MATCHfail;
      return Visitor::VERRORABORT;
    }
    // Key (PRI for primary, MUL if first in a multiple-field key
    // Save primary key info, if any
    if (std::string(colDescrip[3]) == std::string("PRI")) {
      m_primColName = myName;
    }

    // Field 4 is default
    // Extra (may say auto_increment)
    bool autoInc = (colDescrip[5] == "auto_increment");
    if (autoInc != col->isAutoIncrement()) {
      m_matchReturn = MATCHfail;
      return Visitor::VERRORABORT;
    }
    return Visitor::VCONTINUE;
  }

  bool MysqlConnection::checkDType(Datatype* dtype, 
                                   const std::string& sqlType) {
    std::string base;
    int sqlSize;
    if (dtype->getType() != Datatype::TYPEenum) {
      sqlSize = extractSize(sqlType);
    }

    // Cases  char, varchar, enum and datetime are handled entirely within
    // the switch statement, but most do the bulk of the work in
    // common, after the switch.
    switch (dtype->getType()) {
    case Datatype::TYPEenum: {
      base = "enum";
      if (sqlType.find(base) != 0) {
        m_matchReturn = MATCHfail;
        return false;
      }
      Enum* ourEnum = dtype->getEnum();
      // Finally compare local list of choices to those listed in sqlType
      // Local list is a vector; in sqlType they're quoted, comma separated
      return compareEnumList(ourEnum->getChoices(), sqlType);
    }
    case Datatype::TYPEvarchar: {
      base = "varchar";
      if (sqlType.find(base) != 0) {
        m_matchReturn = MATCHfail;
        return false;
      }
      // size in db must be at least as large as size in Col.
      if (sqlSize < dtype->getOutputSize()) {
        m_matchReturn = MATCHfail;
        return false;
      }
      else if (sqlSize > dtype->getOutputSize()) {
        m_matchReturn = MATCHcompatible;
      }
      return true;
    }
    case Datatype::TYPEchar: {
      base = "char";
      if (sqlType.find(base) != 0) {
        m_matchReturn = MATCHfail;
        return Visitor::VERRORABORT;
      }
      //  For char datatype unspecified size is equivalent to size=1
      if (!sqlSize) sqlSize = 1;
      // size in db must be at least as large as size in Col.
      if (sqlSize < dtype->getOutputSize()) {
        m_matchReturn = MATCHfail;
        return false;
      }
      else if (sqlSize > dtype->getOutputSize()) {
        m_matchReturn = MATCHcompatible;
      }
      return true;
    }
    case Datatype::TYPEdatetime: {
      if (sqlType != "datetime") {
        m_matchReturn = MATCHfail;
        return false;
      }
      return true;
    }


    case Datatype::TYPEtimestamp: {
      base = "timestamp";
      break;
    }
    case Datatype::TYPEint: {
      base = "int";
      break;
    }
    case Datatype::TYPEmediumint: {
      base = "mediumint";
      break;
    }
    case Datatype::TYPEsmallint: {
      base = "smallint";
      break;
    }
    case Datatype::TYPEreal: 
    case Datatype::TYPEdouble: {
      base = "double";
      break;
    }
    default: {  // Indicates bad xml file input.  Applications
                //should have exited already
      m_matchReturn = MATCHfail;
      return false;
    }
    }     // end switch
    if (sqlType.find(base) != 0) {
      m_matchReturn = MATCHfail;
      return Visitor::VERRORABORT;
    }
    // Now check size.  It's only for display, so mismatch is not failure
    if (sqlSize != dtype->getOutputSize()) {
      m_matchReturn = MATCHcompatible;
    }

    return true;
  }



  Visitor::VisitorState MysqlConnection::visitIndex(Index* ) {
    return Visitor::VCONTINUE;
    // might put something real here later
  }

  Visitor::VisitorState MysqlConnection::visitAssertion(Assertion*) {
    return Visitor::VCONTINUE;
  }

} // end namespace rdbModel
