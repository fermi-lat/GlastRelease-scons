// $Header$

/*
#ifdef  WIN32
#include <windows.h>
#endif
*/

#include "calibUtil/Metadata.h"
#include "facilities/Util.h"
#include "facilities/Timestamp.h"
#include "rdbModel/Management/Manager.h"
#include "rdbModel/Management/XercesBuilder.h"

#include "rdbModel/Db/MysqlConnection.h"
#include "rdbModel/Db/MysqlResults.h"
#include "rdbModel/Rdb.h"
#include "rdbModel/RdbException.h"
#include "rdbModel/Tables/Assertion.h"
#include "rdbModel/Tables/Table.h"
#include "rdbModel/Tables/Column.h"
#include <iostream>
#include <cstdio>

namespace calibUtil {

  // Might be getting rid of this
  //  const unsigned int Metadata::s_rowReady = Metadata::eOpened 
  //  | Metadata::eValid  | Metadata::eInputDesc | Metadata::eComment;

  Metadata::Metadata(const std::string& host, const std::string& table,
                     const std::string& dbName)  
    : m_readCxt(0), m_writeCxt(0), //  m_row(""), m_rowStatus(0), 
      m_host(host), m_table(table), m_dbName(dbName), m_man(0), m_rdb(0),
      m_match(false) {
    if (table.compare("*") == 0) m_table = std::string("$(MYSQL_METATABLE)");
    if (host.compare("*") == 0) m_host = std::string("$(MYSQL_HOST)");

    int nsub = facilities::Util::expandEnvVar(&m_table);
    // If this doesn't work, use default
    if (nsub < 0) m_table = std::string("metadata");
  }

  Metadata::~Metadata() {
    disconnectRead();
    disconnectWrite();
    if (m_man) delete m_man;
  }


  Metadata::eRet Metadata::fetchUser(std::string& user) {
    // WARNING: Windows and Linux/Solaris use different standard variables
    // for login name.  The test used below isn't quite right
    // since one could conceivably compile with gcc on cygwin/Windows.
#ifdef    __GNUG__ 
    user = std::string("$(USER)");
#else
    user = std::string("$(USERNAME)");
#endif

    int nsub = facilities::Util::expandEnvVar(&user);
    if (nsub == 1) {
      return RETOk;
    }
    else {
      user = std::string("");
      return RETBadValue;
    }
  }

  // The next 5 methods concern connection to the server
  bool Metadata::connect(rdbModel::Connection* cxt, std::string& host,
                         const std::string& user, 
                         const std::string& pw, eRet& err,
                         const std::string& dbName)  {

    int nSub = facilities::Util::expandEnvVar(&host);
    if (nSub < 0) {
      err = RETBadHost;
      return false;
    }

    bool connected = cxt->open(host, user, pw, dbName);
    if (connected) {
      err = RETOk;
      return true;
    }
    else {
      err = RETNoConnect;
      return false;
    }
  }

  bool Metadata::connectRead(eRet& err) {
    if (m_readCxt == 0) {
      // for now use std::cout, std::cerr 
      m_readCxt = new rdbModel::MysqlConnection();
      bool ok = connect(m_readCxt, m_host, std::string("glastreader"), 
                        std::string("glastreader"),  err, m_dbName);
      if (!ok) {
        delete m_readCxt;
        m_readCxt = 0;
      }  else { // look for compatible schema
        std::string schema = 
          std::string("$(RDBMODELROOT)/xml/")+ m_dbName + ".xml"; 
        err = compareSchema(m_readCxt, schema);
      }
      return ok;
    }
    else return true;
  }

  bool Metadata::connectWrite(eRet& err) {
    if (m_writeCxt == 0) {
      m_writeCxt = new rdbModel::MysqlConnection();
      bool ok;
      if (m_dbName == std::string("calib_test") ) {
        ok = connect(m_writeCxt, m_host, std::string("calib_tester"), 
                     std::string("udine"), err, m_dbName);
      }
      else { 
        ok = connect(m_writeCxt, m_host, std::string("calibrator"), 
                     std::string("calibr8tor"), err, m_dbName);
      }
      if (!ok) {
        delete m_readCxt;
        m_readCxt = 0;
      }     else { // look for compatible schema
        std::string schema = 
          std::string("$(RDBMODELROOT)/xml/")+ m_dbName + ".xml"; 
        err = compareSchema(m_writeCxt, schema);
      }

      return ok;
    }
    else return true;
  }

  void Metadata::disconnectRead() {
    if (m_readCxt) {
      m_readCxt->close();
      delete m_readCxt;
      m_readCxt = 0;
    }
  }

  void Metadata::disconnectWrite() {
    if (m_writeCxt) {
      m_writeCxt->close();
      delete m_writeCxt;
      m_writeCxt = 0;
    }
  }

  Metadata::eRet 
  Metadata::findSoonAfter(unsigned int *ser,
                          const std::string& calibType, 
                          const facilities::Timestamp* update_start,
                          const facilities::Timestamp* update_end,
                          unsigned int levelMask, 
                          const std::string& instrument,
                          const std::string& flavor) {
    using namespace rdbModel;

    eRet ret;
    *ser = 0;
    if (!m_readCxt) {
      connectRead(ret);
      if (ret != RETOk) return ret;
    }

    StringVector orderBy;
    orderBy.push_back("update_time");
    

    std::vector<Assertion::Operator *> conditions;
    conditions.reserve(8);

    Assertion::Operator completeOp(OPTYPEequal, "completion", "OK",
                                   FIELDTYPEold, FIELDTYPElit);
    //                               false, true);
    Assertion::Operator instOp(OPTYPEequal, "instrument", instrument,
                               FIELDTYPEold, FIELDTYPElit);
    //                           false, true);
    Assertion::Operator calibTypeOp(OPTYPEequal, "calib_type", calibType,
                                    FIELDTYPEold, FIELDTYPElit);
    //                                false, true);
    Assertion::Operator flavorOp(OPTYPEequal, "flavor", flavor, 
                                 FIELDTYPEold, FIELDTYPElit);
                                 //false, true);
    Assertion::Operator startOp(OPTYPElessOrEqual, update_start->getString(),
                                "update_time", 
                                FIELDTYPEold, FIELDTYPElit);
                                // true, false);
    
    conditions.push_back(&completeOp);
    conditions.push_back(&instOp);
    conditions.push_back(&calibTypeOp);
    conditions.push_back(&flavorOp);
    conditions.push_back(&startOp);

    // Have to handle this one differently since it's optional
    // and Assertion::Operator doesn't have 'set' methods
    Assertion::Operator* endOp = 0;
    if (update_end) {
      endOp = new Assertion::Operator(OPTYPEgreaterOrEqual, 
                                      update_end->getString(),
                                      "update_time", 
                                      FIELDTYPElit, FIELDTYPEold);
                                      //true, false);
      conditions.push_back(endOp);
    }

    ret = doSelect(levelMask, conditions, orderBy, ser);
    if (endOp) delete endOp;
    return ret;
  }

  Metadata::eRet Metadata::findBest(unsigned int *ser,
                          const std::string& calibType, 
                          const facilities::Timestamp& timestamp, 
                          unsigned int levelMask,    // could have default
                          const std::string& instrument, // could have default
                          const std::string& flavor)
  {
    using namespace rdbModel;

    eRet ret;
    *ser = 0;
    if (!m_readCxt) {
      connectRead(ret);
      if (ret != RETOk) return ret;
    }

    // Sort rows by timestamp. Would like most recent first
    StringVector orderBy;
    orderBy.push_back("update_time desc");

    std::vector<Assertion::Operator *> conditions;
    conditions.reserve(8);
    Assertion::Operator completeOp(OPTYPEequal, "completion", "OK",
                                   FIELDTYPEold, FIELDTYPElit);
    //                               false, true);
    Assertion::Operator instOp(OPTYPEequal, "instrument", instrument,
                               FIELDTYPEold, FIELDTYPElit);
    //                           false, true);
    Assertion::Operator calibTypeOp(OPTYPEequal, "calib_type", calibType,
                                    FIELDTYPEold, FIELDTYPElit);
    //                                false, true);
    Assertion::Operator flavorOp(OPTYPEequal, "flavor", flavor, 
                                 FIELDTYPEold, FIELDTYPElit);
                                 //  false, true);
    Assertion::Operator vstartOp(OPTYPEgreaterThan, timestamp.getString(), 
                                 "vstart",
                                 FIELDTYPElit, FIELDTYPEold);
    //                      true, false);

    Assertion::Operator vendOp(OPTYPElessThan, timestamp.getString(), "vend",
                               FIELDTYPElit, FIELDTYPEold);
    //                      true, false);

    conditions.push_back(&completeOp);
    conditions.push_back(&instOp);
    conditions.push_back(&calibTypeOp);
    conditions.push_back(&flavorOp);
    conditions.push_back(&vstartOp);
    conditions.push_back(&vendOp);

    // Finally, set a value for PROC_LEVEL and make the query (ies)
    return doSelect(levelMask, conditions, orderBy, ser);
  }


  Metadata::eRet Metadata::getReadInfo(unsigned int serialNo, 
                                       std::string& dataFmt, 
                                       std::string& fmtVersion,
                                       std::string& filename) {
    using namespace rdbModel;
    eRet ret;
    if (!m_readCxt) {
      connectRead(ret);
      if (ret != RETOk) {
        return ret;
      }
    }
    StringVector colNames(3);
    colNames[0] = "data_fmt";
    colNames[1] = "fmt_version";
    colNames[2] = "data_ident";

    std::string serNoVal;
    facilities::Util::itoa(serialNo, serNoVal);
    Assertion::Operator* serOp = 
      new Assertion::Operator(OPTYPEequal, "ser_no", serNoVal, 
                              FIELDTYPEold, FIELDTYPElit);
                              //false, true);
    //    Assertion whereClause(Assertion::WHENwhere, serOp);
    Assertion whereClause(serOp);
    StringVector orderBy;
    orderBy.clear();
    
    // make the query
    ResultHandle* results = 0;
    try {
      results = m_readCxt->select(m_table, colNames, orderBy, &whereClause);
    }
    catch (RdbException ex) {
      std::cout << ex.getMsg();
      return RETMySQLError;
    }
    
    if (!results) return RETMySQLError;
    if (!results->getNRows() ) return RETBadValue;

    std::vector<std::string> fields;
    results->getRow(fields);
    dataFmt = fields[0];
    fmtVersion = fields[1];
    filename = fields[2];

    if ((dataFmt == "") || (filename == "")) return RETBadValue;
    return RETOk;
  }

Metadata::eRet Metadata::getInterval(unsigned int serialNo, 
                                       facilities::Timestamp*& since,
                                       facilities::Timestamp*& till) {
    using namespace rdbModel;
    StringVector orderBy;
    orderBy.clear();

    eRet ret;
    if (!m_readCxt) {
      connectRead(ret);
      if (ret != RETOk) return ret;
    }

    StringVector colNames(2);
    colNames[0] = "vstart";
    colNames[1] = "vend";

    std::string serNoVal;
    facilities::Util::itoa(serialNo, serNoVal);
    Assertion::Operator* serOp = 
      new Assertion::Operator(OPTYPEequal, "ser_no", serNoVal, 
                              FIELDTYPEold, FIELDTYPElit);
                              //false, true);

    // Assertion::Operator serOp(OPTYPEequal, "ser_no", serNoVal, false, true);
    Assertion whereClause(serOp);
    StringVector orderCols;
    orderCols.clear();
    
    ResultHandle* results = 0;

    // make the query
    try {
      results = m_readCxt->select(m_table, colNames, orderBy, &whereClause);
    }
    catch (RdbException ex) {
      std::cout << ex.getMsg();
      // return heap memory
      return RETMySQLError;
    }
    
    if (!results) return RETMySQLError;
    if (!results->getNRows() ) return RETBadValue;

    std::vector<std::string> fields;

    results->getRow(fields);

    since = new facilities::Timestamp(fields[0]);
    till =  new facilities::Timestamp(fields[1]);
    return RETOk;
  }

  /*

| Field       | Type               | Null | Key | Default | How set          |
+-------------+--------------------+------+-----+---------+------------------+
| ser_no      | mediumint(9)       |      | PRI | NULL    | auto_increment   |
| instrument  | varchar(16)        |      |     |         |   [openRecord]   |
| calib_type  | varchar(20)        |      |     |         |   [openRecord    |
| data_fmt    | varchar(10)        |      |     |         |   [openRecord]   |
| data_size   | int(11)            | YES  |     | NULL    |   [optional]     |
| vstart      | datetime           | YES  |     | NULL    |[addValidInterval]|
| vend        | datetime           | YES  |     | NULL    |[addValidInterval]|
| enter_time  | timestamp(14)      | YES  |     | NULL    | automatic        |
| fmt_version | varchar(12)        | YES  |     | NULL    |   [openRecord]   |
| completion  | enum('OK','INC','ABORT')           
                                   | YES  | MUL | NULL    |   [openRecord]   |
| proc_level  | enum('PROD','TEST','DEV', 'SUPSED')
                                   |      |     | TEST    |   [openRecord]   |
| creator     | varchar(255)       | YES  |     | NULL    |   [addCreator]   |
| uid         | varchar(12)        |      |     |         |[insertRecord/
                                                             addUser]        |
| data_ident  | varchar(255)       |      |     |         |   [openRecord]   |
| input_desc  | varchar(255)       | YES  |     | NULL    |[addInputDesc]    |
| notes       | varchar(255)       | YES  |     | NULL    |   [addNotes]     |
+-------------+--------------------+------+-----+---------+------------------+

  */
  /*               Private utilities                           */

  bool Metadata::fetchLevel(std::string& val, unsigned int *levelMask) {
    bool ret = true;
    unsigned int testBit = 0;
    if ((testBit = LEVELProd) & *levelMask) {
      val = std::string("PROD");
    }
    else if ((testBit = LEVELDev) & *levelMask) {
      val = std::string("DEV");
    }
    else if ((testBit = LEVELTest) & *levelMask) {
      val = std::string("TEST");
    }
    else if ((testBit = LEVELSuperseded) & *levelMask) {
      val = std::string("SUPSED");
    }
    else {      // All that's left are unknown bits.  
      testBit = *levelMask;
      std::cerr << "Improper bit(s) in level mask " << testBit 
                << std::endl;
      ret = false;
    }
    *levelMask &= (~testBit);
    return ret;
  }    



  Metadata::eRet 
  Metadata::doSelect(unsigned levelMask, 
                     std::vector<rdbModel::Assertion::Operator *>& conditions, 
                     rdbModel::StringVector& orderBy,
                     unsigned*& ser) {
    using namespace rdbModel;
    *ser = 0;

    StringVector colNames(1);
    colNames[0] = "ser_no"; 

    while (levelMask) {
      std::string val;
      Assertion::Operator* levelOp = 0;
      bool ok = fetchLevel(val, &levelMask);
      if (!ok) {
        std::cerr << "In Metadata::doSelect(..) bad levelMask arg" 
                  << std::endl;
        return RETBadValue;
      }
      levelOp = 
        new Assertion::Operator(OPTYPEequal, val, "PROC_LEVEL", 
                                FIELDTYPElit, FIELDTYPEold);
                                //true, false);
      conditions.push_back(levelOp); 
      
      // make the combined operator
      Assertion::Operator andOp(OPTYPEand, conditions, true);

      //Following creates an assertion such that creator (us) continues
      // to own the associated operator.
      Assertion whereClause(&andOp, 0, true);
      //      Assertion whereClause(Assertion::WHENwhere, &andOp, 0, true);
        
      ResultHandle* results = 0;

      try {               // make the query
        results = m_readCxt->select(m_table, colNames, orderBy, &whereClause);
      }
      catch (RdbException ex) {
        std::cout << ex.getMsg();
        // return heap memory
        delete levelOp;  levelOp = 0;
      }

      if (!results) {   // Error. Should have ResultHandle even if 0 rows.
        return RETMySQLError;
      }

      // handle result.  
      if (results->getNRows() == 0) {
        conditions.pop_back(); 
        if (levelOp) {
          delete levelOp;
          levelOp = 0;
        }
        continue;   // try next proc level
      }
      std::vector<std::string> fields;
      results->getRow(fields);
      *ser = facilities::Util::stringToInt(fields[0]);
      return RETOk;
    }
    return RETOk;
  }

  int  Metadata::registerCalib(const std::string& inst, 
                               const std::string& flavor,
                               const std::string& calib_type, 
                               const std::string& data_ident,
                               const std::string& data_fmt,
                               const facilities::Timestamp& vstart,
                               const facilities::Timestamp& vend,
                               const std::string& input_desc, 
                               const std::string& notes,
                               const std::string& proc_level,
                               const std::string& locale,
                               const std::string& fmt_version,
                               const std::string& completion,
                               const std::string& input_start,
                               const std::string& input_end) {
    using namespace rdbModel;

    eRet ret;

    if (vend < vstart) {
      std::cout << 
        "Metadata::registerCalib: Error in metadata. vend < vstart " <<
        std::endl;
      return 0;
    }
    if (!m_writeCxt) {
      connectWrite(ret);
      if (ret != RETOk) return 0; // we or connectWrite should throw exception
    }
    StringVector cols;
    StringVector vals;
    StringVector nullCols;

    cols.reserve(24);
    vals.reserve(24);
    nullCols.reserve(16);

    if (inst.size() * calib_type.size() * flavor.size() * data_fmt.size() 
        * data_ident.size() * proc_level.size() * completion.size() 
        * locale.size() == 0) {   // something is null that shouldn't be
      return 0;      // should perhaps throw exception
    }
    cols.push_back("instrument"); vals.push_back(inst);
    cols.push_back("calib_type"); vals.push_back(calib_type);
    cols.push_back("flavor"); vals.push_back(flavor);
    cols.push_back("data_fmt"); vals.push_back(data_fmt);
    cols.push_back("data_ident"); vals.push_back(data_ident);
    cols.push_back("vstart"); vals.push_back(vstart.getString());
    cols.push_back("vend"); vals.push_back(vend.getString());
    cols.push_back("proc_level"); vals.push_back(proc_level);
    cols.push_back("completion"); vals.push_back(completion);
    cols.push_back("locale"); vals.push_back(locale);

    //  These, however, may be null
    if (input_desc.size() > 0 ) {
      cols.push_back("input_desc"); vals.push_back(input_desc); 
    }  else nullCols.push_back("input_desc");

    if (notes.size() > 0 ) {
      cols.push_back("notes"); vals.push_back(notes);
    }    else nullCols.push_back("notes");

    if (fmt_version.size() > 0 )
    {
      cols.push_back("fmt_version"); 
      vals.push_back(fmt_version);
    } 
    //    else nullCols.push_back("fmt_version");

    if (input_start.size() > 0) {
      cols.push_back("input_start"); vals.push_back(input_start);
    }    else nullCols.push_back("input_start");
    
    if (input_end.size() > 0) {
      cols.push_back("input_end"); vals.push_back(input_end);
    }    else nullCols.push_back("input_end");

    // The service -- that's us -- is responsible for creator, uid, enter_time
    cols.push_back("creator"); vals.push_back("Metadata::registerCalib");
    std::string uid;
    fetchUser(uid);
    cols.push_back("uid"); vals.push_back(uid);
    facilities::Timestamp curTime;
    cols.push_back("enter_time");vals.push_back(curTime.getString());
    // update_time is set automatically by MySQL, but MySQL uses
    // local timezone rather than gmt, so we have set it explicitly
    cols.push_back("update_time");vals.push_back(curTime.getString());


    if (m_rdb) {
      bool ok = checkValues(cols, vals);
      if (ok) checkNulls(nullCols);
      if (!ok) return 0;
    }

    // ser_no gets set automatically by MySQL
    int ser_no;
    if (!(m_writeCxt->insertRow(m_table, cols, vals, &ser_no, &nullCols)) ) {
      return 0;
    }    else {
      adjustVend(ser_no);
      return ser_no;
    }
  }

  Metadata::eRet Metadata::compareSchema(rdbModel::Connection* conn,
                                         const std::string& schema) {
    using namespace rdbModel;


    if (m_man) { // already did this
      return (m_match) ? RETOk : RETNoSchemaMatch;
    }
    m_man = rdbModel::Manager::getManager();

    m_man->setBuilder(new rdbModel::XercesBuilder);
    m_man->setInputSource(schema);

    // Maybe first check if file exists?  E.g., try opening for read

    // good errcode is 0
    int errcode = m_man->build();

    if (errcode) {
      std::cerr << "Error in database description file " << schema 
                << std::endl;
      std::cerr << "Parse failed with error " << errcode << std::endl;
      return RETBadCnfFile;
    }
    m_rdb = m_man->getRdb();

    rdbModel::MATCH match = conn->matchSchema(m_rdb, false);

    switch (match) {
    case rdbModel::MATCHequivalent:
    case rdbModel::MATCHcompatible:
      m_match = true;
      return RETOk;
    case rdbModel::MATCHfail:
      std::cout << "XML schema and MySQL database are NOT compatible" 
                << std::endl;
      return RETBadCnfFile;
    case rdbModel::MATCHnoConnection:
      std::cout << "Connection failed while attempting match" << std::endl;
      return RETNoConnect;
    }
    return RETBadValue;
  }

  bool Metadata::checkValues(const rdbModel::StringVector& cols,
                             const rdbModel::StringVector& vals) const {
    unsigned nCol = cols.size();
    rdbModel::Table* table = m_rdb->getTable(m_table);

    
    for (unsigned iCol = 0; iCol < nCol; iCol++) {
      rdbModel::Column* col = table->getColumnByName(cols[iCol]);
      if (!col->okValue(vals[iCol])) {
        std::cerr << "Value " << vals[iCol] << " not allowed for column "
                  << cols[iCol] << " in table " << m_table << std::endl;
        return false;
      }
    }
    return true;
  }

  bool Metadata::checkNulls(const rdbModel::StringVector& cols) const {
    unsigned nCol = cols.size();
    rdbModel::Table* table = m_rdb->getTable(m_table);

    
    for (unsigned iCol = 0; iCol < nCol; iCol++) {
      rdbModel::Column* col = table->getColumnByName(cols[iCol]);
      if (!col->nullAllowed()) {
        std::cerr << "Column "
                  << cols[iCol] << " in table " << m_table
                  << " is not nullable" << std::endl;
        return false;
      }
    }
    return true;
  }


  unsigned Metadata::adjustVend(int newSer) {
    using namespace rdbModel;

    StringVector getCols;
    StringVector orderBy;

    orderBy.clear();
    getCols.reserve(7);

    std::string serString;
    facilities::Util::itoa(newSer, serString);

    getCols.push_back("flavor");
    getCols.push_back("calib_type");
    getCols.push_back("completion");
    getCols.push_back("instrument");
    getCols.push_back("proc_level");
    getCols.push_back("vstart");

    ResultHandle* results = 0;
    eRet err;
    Assertion* where = 0;
    try {
      if (!m_writeCxt) {
        if (!connectWrite(err)) return 0;
      }

      Assertion::Operator* serOp = 
        new Assertion::Operator(OPTYPEequal, "ser_no", serString, 
                                FIELDTYPEold, FIELDTYPElit);
                                //false, true);
    
      where = new Assertion(serOp);

      // Fetch information for new row: vstart, flavor, completion,
      // proc_level, calib_type, flavor
      results = m_writeCxt->select(m_table, getCols, orderBy, where);
      delete where;
      where = 0;
    }
    catch (RdbException ex) {
      std::cout << ex.getMsg();
      delete where;         // return heap memory
      return 0;
    }
    if (!results) { // This is an error. Should be non-null even if no rows
      std::cout << "MySQL failure in SELECT" << std::endl;
      return 0;  // nothing to fix
    }
    if (results->getNRows() != 1) { // also a problem
      std::cout << "Look-up of serial# " << serString << " failed"
                << std::endl;
      return 0;
    }
    std::vector<std::string> fields;
    results->getRow(fields);
    if (fields[2] != "OK") return 0;  // don't bother fixing in this case

    // Now do an update on rows satisfying
    // ((flavor="f") && (calib_type = "c") && (completion= "OK") && 
    // (instrument = "i") &&  (proc_level = "p") && 
    //  (vstart < "new-start" (vend > "new-vstart") );
    std::vector<Assertion::Operator *> conditions;
    conditions.reserve(7);
    for (unsigned ix = 0; ix < 5; ix++) {
      conditions.push_back(new Assertion::Operator(OPTYPEequal, getCols[ix],
                                                   fields[ix], 
                                                   FIELDTYPEold, FIELDTYPElit));
                                                   // false, true));
    }
    conditions.push_back(new Assertion::Operator(OPTYPElessThan, "vstart",
                                                 fields[5], 
                                                 FIELDTYPEold, FIELDTYPElit));
    //false, true));
    conditions.push_back(new Assertion::Operator(OPTYPEgreaterThan, "vend",
                                                 fields[5], 
                                                 FIELDTYPEold, FIELDTYPElit));
                                                 // false, true));
                                            
    Assertion::Operator* andOp = 
      new Assertion::Operator(OPTYPEand, conditions);
    where = new Assertion(andOp);
    
    StringVector toUpdate;
    toUpdate.push_back("vend");
    StringVector newVal;
    newVal.push_back(fields[5]);
    // also update update_time. If we leave it to MySQL, won't be GMT
    facilities::Timestamp curTime;
    toUpdate.push_back("update_time");newVal.push_back(curTime.getString());

    unsigned nModified = m_writeCxt->update(m_table, toUpdate, newVal, where);
    delete where;
    return nModified;
  }

}
