// $Header$


#include "calibUtil/Metadata.h"
#include "facilities/Util.h"
#include "mysql.h"
#include <strstream>
#include <cstdio>

namespace calibUtil {
  // Initialize static members
  // Consider changing all this staticness to avoid chance of difficulty
  // when building calibUtil shareable (for now it's just a static library).
  //  MYSQL* Metadata::readCxt = 0;
  //  MYSQL* Metadata::writeCxt = 0;
  //  bool Metadata::writeConnect = false;
  //  bool Metadata::readConnect = false;
  //  MYSQL     Metadata::readCxt;
  //  MYSQL     Metadata::writeCxt;

  const unsigned int Metadata::rowReady = Metadata::eOpened | Metadata::eValid 
  | Metadata::eInputDesc | Metadata::eComment;


  Metadata::Metadata() : readCxt(0), writeCxt(0), row(""), rowStatus(0), 
    table("$(MYSQL_METATABLE)") {
    int nsub = facilities::Util::expandEnvVar(&table);
    // IF this doesn't work, use default
    if (!nsub) table = std::string("metadata");
  }

  Metadata::~Metadata() {
    disconnectRead();
    disconnectWrite();
  }


  Metadata::eRet Metadata::addCreator(std::string creator) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eCreator) return RETWrongState;

    row += ",creator='"; row += creator; row += "'";
    rowStatus |= eCreator;
    return RETOk;
  }

  //  void Metadata::addField(const std::string& fieldName, 
  //                          const std::string& fieldValue) {
  // To be written
  //                       or maybe not
  //  }

  Metadata::eRet Metadata::addInputDesc(std::string desc) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eInputDesc) return RETWrongState;

    row += ",input_desc='"; row += desc; row+="'";
    rowStatus |= eInputDesc;
    return RETOk;
  }

  Metadata::eRet Metadata::addNotes(std::string comment) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eComment) return RETWrongState;

    row += ",notes='"; row += comment; row += "'";
    rowStatus |= eComment;
    return RETOk;
  }

  // WARNING: Windows and Linux/Solaris use different standard variables
  // for login name.  The test used below isn't quite right
  // since one could conceivably compile with gcc on cygwin/Windows.
  Metadata::eRet Metadata::addUser() {
#ifdef    __GNUG__ 
    std::string user("$(USER)");
#else
    std::string user("$(USERNAME)");
#endif
    row +=",uid='"; 

    int nsub = facilities::Util::expandEnvVar(&user);
    if (nsub == 1) {
      row += user; row+="'";
      return RETOk;
    }
    else {
      row += "'unkown'";
      return RETBadValue;
    }
  }

  Metadata::eRet Metadata::addValidInterval(Timestamp startTime, 
                                            Timestamp endTime) {
    if (!(rowStatus & eOpened) ) return RETWrongState;

    if (rowStatus & eValid) return RETWrongState;
    
    row += ", vstart='"; row += startTime.timeString();
    row += "', vend='"; row += endTime.timeString();
    row += "'";
    rowStatus |=  eValid;
    return RETOk;   // or something else
  }

  // Explicit clear of in-progress record building
  void Metadata::clearRecord() {
    row.resize(0);
    rowStatus = 0;
  }

  // The next 5 methods concern connection to the server
  bool Metadata::connect(MYSQL* cxt, const std::string& user, 
                         const std::string& pw, eRet& err)  {

    std::string host("$(MYSQL_HOST)");

    int nSub = facilities::Util::expandEnvVar(&host);
    if (!nSub) {
      err = RETBadHost;
      return false;
    }

    mysql_init(cxt);
    MYSQL *connected;
    connected = mysql_real_connect(cxt, host.c_str(), user.c_str(),
                                   pw.c_str(), "calib", 0, NULL, 0);
    if (connected != 0) {
      cxt = connected;
      err = RETOk;
      return true;
    }
    else {
      err = RETNoConnect;
      return false;
    }
  }

  bool Metadata::connectRead(eRet& err) {
    if (readCxt == 0) {
      readCxt = new MYSQL;
      return connect(readCxt, std::string("glastreader"), 
                      std::string("glastreader"), err);
    }
    else return true;
  }

  bool Metadata::connectWrite(eRet& err) {
    if (writeCxt == 0) {
      writeCxt = new MYSQL;
      return connect(writeCxt, std::string("calibrator"), 
                      std::string("calibr8tor"), err);
    }
    else return true;
  }


  void Metadata::disconnectRead() {
    if (readCxt) {
      mysql_close(readCxt);
      delete readCxt;
      readCxt = 0;
    }
  }

  void Metadata::disconnectWrite() {
    if (writeCxt) {
      mysql_close(writeCxt);
      delete writeCxt;
      writeCxt = 0;
    }
  }

  Metadata::eRet Metadata::findBest(unsigned int *ser,
                                    eCalibType calibType, 
                                    const Timestamp& timestamp,
                                    unsigned int levelMask, 
                                    eInstrument instrument) {
    const std::string* const ctypeStr = getCalibTypeStr(calibType);
    const std::string* const instStr = getInstrumentStr(instrument);
    if ((ctypeStr != 0) && (instStr !=0)) {
      return findBest(ser, *ctypeStr, timestamp, levelMask,
                      *instStr);
    }
    else return RETBadValue;
  }

  Metadata::eRet Metadata::findBest(unsigned int *ser,
                          const std::string& calibType, 
                          const Timestamp& timestamp, 
                          unsigned int levelMask,    // could have default
                          const std::string& instrument) // could have default
  {
    eRet ret;
    *ser = 0;
    if (!readCxt) {
      connectRead(ret);
      if (ret != RETOk) return ret;
    }

    // Sort rows by timestamp, most recent first
    std::string 
      query("select ser_no,vstart,vend, enter_time,completion from ");
    query += table;
    query += " where ";
    query += "completion = 'OK' and instrument ='"; query += instrument;
    query += "' and calib_type ='";
    query += calibType;
    query += "' and '";
    query += timestamp.timeString();
    query += "'> vstart and vend > '";
    query += timestamp.timeString();
    query += "'";
    while (levelMask) {
      // find highest priority known bit in levelMask, add appropriate
      // clause onto query.  Return false if no known bits found
      std::string q = query;
      bool ok = addLevel(q, &levelMask);  
      if (!ok) {
        std::cerr << "In Metadata::findBest(..) bad levelMask arg" 
                  << std::endl;
        return RETBadValue;
      }
      q += " order by enter_time desc ";

      int ret = mysql_query(readCxt, q.c_str());
      if (ret) {
        std::cerr << "MySQL error during SELECT, code " << ret << std::endl;
        std::cerr << "query was: " << std::endl;
        std::cerr << q << std::endl;
        return RETMySQLError;
      }
      
      MYSQL_RES *myres = mysql_store_result(readCxt);

      // Since we're doing a query, a result set should be returned
      // even if there are no rows in the result.  A null pointer 
      // indicates an error.
      if (!myres) {
        return RETMySQLError;
      }
      else if (mysql_num_rows(myres) ) {       // Get serial # of first row
        MYSQL_ROW myRow = mysql_fetch_row(myres);
        // serial number is pointed to by myRow[0]
        *ser = atoi(myRow[0]);
        return RETOk;
      }
      // otherwise there was no error, but also no rows returned
      // which matched the query, so keep going.
    }


    return RETOk;

  }

  // Not done yet, but this is more than enough for testing.
  const std::string* const Metadata::getCalibTypeStr(eCalibType cType) {
    static std::string acdEff("ACD_Eff");
    static std::string tkrBadChan("TKR_Bad_Chan");
    static std::string calLightAtt("CAL_Light_Att");
    static std::string calLightAsym("CAL_Light_Asym");
    switch(cType) {
    case CTYPE_ACDEff:
      return &acdEff;
    case CTYPE_TKRBadChan:
      return &tkrBadChan;
    case CTYPE_CALLightAtt:
      return &calLightAtt;
    case CTYPE_CALLightAsym:
      return &calLightAsym;
    default:
      return 0;
    }

  }

  const std::string* const Metadata::getDataFmtStr(eDataFmt fmt) {
    static std::string xml("XML");
    static std::string root("ROOT");
    switch (fmt) {
    case FMTXml:
      return &xml;
    case FMTRoot:
      return &root;
    default:
      return 0;
    }
  }

  const std::string* const Metadata::getInstrumentStr(eInstrument inst) {
    static std::string lat("LAT");
    static std::string btem("BTEM");
    static std::string bfem("BFEM");
    static std::string em("EM");
    static std::string fu("FU");
    static std::string cu("CU");
    
    switch(inst) {
    case INSTLat:
      return &lat;
    case INSTBtem:
      return &btem;
    case INSTBfem:
      return &bfem;
    case INSTEm:
      return &em;
    case INSTFu:
      return &fu;
    case INSTCu:
      return &cu;
    default:
      return 0;
    }
  }

  Metadata::eRet Metadata::getReadInfo(unsigned int serialNo, 
                             eDataFmt&    dataFmt, 
                             std::string& fmtVersion,
                             std::string& filename) {
    eRet ret;
    if (!readCxt) {
      connectRead(ret);
      if (ret != RETOk) return ret;
    }
    std::string q("select data_fmt, fmt_version, data_ident from ");
    q += table;
    q += " where ser_no=";

    char serBuf[20];
    sprintf(serBuf, "%i", serialNo);

    q += serBuf;

    int myRet = mysql_query(readCxt, q.c_str());
    if (myRet) {
      std::cerr << "MySQL error during SELECT, code " << ret << std::endl;
      return RETMySQLError;
    }
      
    MYSQL_RES *myres = mysql_store_result(readCxt);
    
    if (mysql_num_rows(myres) ) {  // must have been a good serial number
      MYSQL_ROW myRow = mysql_fetch_row(myres);
      std::string fmtString(myRow[0]);
      
      if (!(fmtString.compare(std::string("XML"))) ) {
        dataFmt = FMTXml;
      }
      else if (!(fmtString.compare(std::string("ROOT"))) ) {
        dataFmt = FMTRoot;
      }
      else dataFmt = FMTUnknown;


      fmtVersion = std::string(myRow[1]);
      filename = std::string(myRow[2]);
      return RETOk;
    }
    else return RETBadValue;
  }

  /*
    openRecord, insertRecord, and various add.. methods are used to
    write a record to the database.  The table below lists all columns
    and indicates how they get their values.

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
  Metadata::eRet Metadata::insertRecord(unsigned int* serialNo) {
    // Send it off (To be written)

    if (serialNo) *serialNo = 0;
    if ((!(rowStatus & rowReady)) == rowReady) return RETWrongState;

    if (!(rowStatus & eCreator)) {
      addCreator("calibUtil::Metadata::insertRecord");
    }
    addUser();
    // Actually write it...
    int ret = mysql_query(writeCxt, row.c_str());
    clearRecord();

    if (ret) {
      std::cerr << "MySQL error during INSERT, code " << ret << std::endl;
      return RETMySQLError;
    }
    if (serialNo) {
      *serialNo = mysql_insert_id(writeCxt);
    }
    return RETOk;   // or something else if we failed 
  }


  Metadata::eRet Metadata::openRecord(eInstrument        inst,
                                      eCalibType         ctype,
                                      eDataFmt           fmt,
                                      const std::string& fmtVersion,
                                      const std::string& dataIdent, 
                                      eCompletion        completion,
                                      eLevel             procLevel){
    if (row.size() > 0) {
      std::cerr << 
        "calibration metadata record build already in progress" << std::endl;
      return RETWrongState;
    }
    // Do checking in advance of column values which come from an
    // enumerated sets.
    const std::string* const  instString = getInstrumentStr(inst);
    const std::string* const  calibTypeStr = getCalibTypeStr(ctype);
    const std::string* const  fmtString = getDataFmtStr(fmt);
    const std::string* const  cmpString  = getCompletionStr(completion);
    const std::string* const  levelString = getProcLevelStr(procLevel);

    if (!fmtString || !instString || !calibTypeStr ||
        !cmpString || !levelString) return RETBadValue;

    return openRecord(*instString, *calibTypeStr, *fmtString, fmtVersion,
                      dataIdent, *cmpString, *levelString);
  }

  /// Similar to above but with (nearly) all string arguments for
  /// maximum flexibility, less argument validation
  Metadata::eRet Metadata::openRecord(const std::string& inst, 
                                      const std::string& calibType,
                                      const std::string& dataFmt,
                                      const std::string& fmtVersion,
                                      const std::string& dataIdent, 
                                      const std::string& completion,
                                      const std::string& level) {

    if (row.size() > 0) {
      std::cerr << 
        "calibration metadata record build already in progress" << std::endl;
      return RETWrongState;
    }

    if (!writeCxt) {
      eRet ret;
      connectWrite(ret);
      if (ret != RETOk) return ret;
    }

    //    += is a synonym for row.append(..)
    row += "insert into "; row += table;
    row += " set instrument='";    row += inst;
    row += "', calib_type='";  row += calibType;
    row += "', data_fmt='"; row += dataFmt;
    
    row += "', fmt_version='"; row += fmtVersion;
    row += "', completion='"; row += completion;
    row += "', data_ident='"; row += dataIdent;
    row += "', proc_level='"; row += level; row+="'";
    rowStatus = eOpened;
    return RETOk;
    }
  /*               Private utilities                           */

  // Internal utility. Sticks on the last part of the "where" clause
  // -- the part for dev/prod/..etc.
  // Return "true" if we added something to the query (found a bit
  // we understand in *levelMask) otherwise return false.
  bool Metadata::addLevel(std::string& q, unsigned int *levelMask) {
    
    bool ret = true;
    unsigned int testBit = 0;
    if ((testBit = LEVELProd) & *levelMask) {
      q += "and (proc_level='PROD')";
    }
    else if ((testBit = LEVELDev) & *levelMask) {
      q += "and (proc_level='DEV')";
    }
    else if ((testBit = LEVELTest) & *levelMask) {
      q += "and (proc_level='TEST')";
    }
    else if ((testBit = LEVELSuperseded) & *levelMask) {
      q += "and (proc_level='SUPSED')";
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

  // Internal utilities to do validation of certain fields.  On the
  // whole it's better to do this ourselves rather than depending on
  // MySQL to do it.  The latter is possible in some cases, but 
  // handling and reporting errors is difficult.

  const std::string* const Metadata::getCompletionStr(eCompletion cmpl) {
    static std::string ok("OK");
    static std::string inc("INC");
    static std::string abort("ABORT");
    switch (cmpl) {
    case CMPLOk:
      return &ok;
    case CMPLInc:
      return &inc;
    case CMPLAbort:
      return &abort;
    default:
      return 0;
    }
  }

  const std::string* const Metadata::getProcLevelStr(eLevel level) {
    static std::string prod("PROD");
    static std::string dev("DEV");
    static std::string test("TEST");
    static std::string sup("SUPSED");
    switch(level) {
    case LEVELProd:
      return &prod;
    case LEVELDev:
      return &dev;
    case LEVELTest:
      return &test;
    case LEVELSuperseded:
      return &sup;
    default:
      return 0;
    }
  }

}
