// $Header$


#include "calibUtil/Metadata.h"
#include "facilities/Util.h"
#include "mysql.h"
#include <strstream>


namespace calibUtil {
  // Initialize static members
  // Consider changing all this staticness to avoid chance of difficulty
  // when building calibUtil shareable (for now it's just a static library).
  MYSQL* Metadata::readCxt = 0;
  MYSQL* Metadata::writeCxt = 0;

  const unsigned int Metadata::rowReady = Metadata::eOpened | Metadata::eValid 
  | Metadata::eInputDesc | Metadata::eComment;


  Metadata::Metadata() : row(""), rowStatus(0), table("$(MYSQL_METATABLE)") {
    int nsub = facilities::Util::expandEnvVar(&table);
    // IF this doesn't work, use default
    if (!nsub) table = std::string("metadata");
  }

  // The next 5 methods concern connection to the server
  bool Metadata::connectRead(eRet& err) {
    return connect(readCxt, std::string("reader"), err);
  }

  bool Metadata::connectWrite(eRet& err) {
    return connect(writeCxt, std::string("calibrator"), err);
  }


  bool Metadata::connect(MYSQL* cxt, const std::string& grp, eRet& err)  {
    std::string cnfFile("$(CALIBUTILROOT)/src/.mysqlopt.cnf");

    // Translate environment variables
    int nSub = facilities::Util::expandEnvVar(&cnfFile);

    if (nSub != 1) {
      err = RETBadCnfFile;
      return false;
    }

    std::string host("$(MYSQL_HOST)");

    nSub = facilities::Util::expandEnvVar(&host);
    if (!nSub) {
      err = RETBadHost;
      return false;
    }


    mysql_init(cxt);

    mysql_options(cxt, MYSQL_READ_DEFAULT_FILE, cnfFile.c_str());
    mysql_options(cxt, MYSQL_READ_DEFAULT_GROUP, grp.c_str());

    MYSQL *connected = mysql_real_connect(cxt, host.c_str(), 0, 0, "calib",
                                 0, (char *) NULL, 0);
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

  void Metadata::disconnectRead() {
    mysql_close(readCxt);
    readCxt = 0;
  }

  void Metadata::disconnectWrite() {
    mysql_close(writeCxt);
    writeCxt = 0;
  }


  /*
    The following collection of methods all have to do with writing
    a record to the database.  The table below lists all columns
    and indicates how they get their values.

| Field               | Type           | Null | Key | Default | How set       |
+---------------------+----------------+------+-----+---------+---------------+
| serial_num          | int(11)        |      | PRI | NULL    | auto_increment|
| instrument          | varchar(20)    | YES  |     | NULL    |  [open]       |
| calib_type          | varchar(40)    | YES  |     | NULL    |  [open]       |
| calib_type_num      | int(11)        | YES  |     | NULL    |  [unused]     |
| data_format         | varchar(20)    |      |     |         |  [open]       |
| format_version      | varchar(20)    | YES  |     | NULL    |  [open]       |
| data_id             | varchar(40)    | YES  |     | NULL    |  [optional]   |
| data_size           | int(11)        | YES  |     | NULL    |  [optional]   |
| validity_start_time | datetime       | YES  |     | NULL    | [addValid]    |
| validity_end_time   | datetime       | YES  |     | NULL    | [addValid]    |
| proc_time           | datetime       | YES  |     | NULL    | [auto]        |
| proc_level          | enum('TEST','DEVELOPMENT','PRODUCTION','SUPERSEDED')
                                       | YES  |     | NULL    |  [open]       |
| calib_status        | enum('OK','INCOMPLETE','ABORTED')        [open]
                                       | YES  |     | NULL    |               |
| creator             | varchar(40)    | YES  |     | NULL    | [defaults]    |
| input_desc          | varchar(255)   | YES  |     | NULL    | [addInDesc]   |
| calib_file_name     | varchar(255)   |      |     |         |  [open]       |
| comments            | varchar(255)   | YES  |     | NULL    |  [addComment] |
+---------------------+----------------+------+-----+---------+---------------+

  */

  Metadata::eRet Metadata::openRecord(const std::string& instr, 
                                      const std::string& calibType,
                                      const std::string& dataFmt, 
                                      const std::string& fmtVersion,
                                      const std::string& dataIdent, 
                                      const std::string& completion,
                                      const std::string& procLevel){
    eRet ret;

    if (!writeCxt) {
      connectWrite(ret);
      if (ret != RETOk) return ret;
    }

    if (row.size() > 0) {
      std::cerr << 
        "calibration metadata record build already in progress" << std::endl;
      return RETWrongState;
    }
    // Do checking in advance of column values which come from an
    // enumerated sets.
    if (!(checkCompletionInput(completion) || 
          checkProcLevelInput(procLevel)) )
      return RETBadValue;

    //    += is a synonym for row.append(..)
    row += "insert into "; row += table;
    row += " set instrument='";    row += instr;
    row += "', calib_type='";  row += calibType;
    row += "', data_fmt='"; row += dataFmt;
    row += "', fmt_version='"; row += fmtVersion;
    row += "', completion='"; row += completion;
    row += "', data_ident='"; row += dataIdent;
    row += "', proc_level='"; row += procLevel; row+="'";
    rowStatus = eOpened;
    return RETOk;
  }



  Metadata::eRet Metadata::insertRecord() {
    // Send it off (To be written)

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
    return RETOk;   // or something else if we failed 
  }

  // Explicit clear of in-progress record building
  void Metadata::clearRecord() {
    row.resize(0);
    rowStatus = 0;
  }

  //  void Metadata::addField(const std::string& fieldName, 
  //                          const std::string& fieldValue) {
  // To be written
  //                       or maybe not
  //  }

  Metadata::eRet Metadata::addValidInterval(Timestamp startTime, 
                                            Timestamp endTime) {
    if (!(rowStatus & eOpened) ) return RETWrongState;

    if (rowStatus & eValid) return RETWrongState;
    
    row += ", validity_start_time='"; row += startTime.timeString();
    row += "', validity_end_time='"; row += endTime.timeString();
    row += "'";
    rowStatus |=  eValid;
    return RETOk;   // or something else
  }

  Metadata::eRet Metadata::addNotes(std::string comment) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eComment) return RETWrongState;

    row += ",notes='"; row += comment; row += "'";
    rowStatus |= eComment;
    return RETOk;
  }

  Metadata::eRet Metadata::addInputDesc(std::string desc) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eInputDesc) return RETWrongState;

    row += ",input_desc='"; row += desc; row+="'";
    rowStatus |= eInputDesc;
    return RETOk;
  }

  Metadata::eRet Metadata::addCreator(std::string creator) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eCreator) return RETWrongState;

    row += ",creator='"; row += creator; row += "'";
    rowStatus |= eCreator;
    return RETOk;
  }

  // WARNING: Need to look into what's available on Windows.
  // Might have to do something different to get user (or just punt)
  Metadata::eRet Metadata::addUser() {
    std::string user("$(USER)");
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

  Metadata::eRet Metadata::findBest(unsigned int *ser,
                          const std::string& calibType, 
                          Timestamp timestamp, 
                          unsigned int levelMask,    // could have default
                          const std::string& instrument) // could have default
  {
    // To be written
    // Retrieve all rows 
    //     * of correct calibType
    //     * for which timestamp is within valid interval
    //     * which have a match with levelMask
    //  Don't need all columns of the rows.  Will need at least
    //       - serial number
    //       - valid start and end times
    //       - proc_time
    //       - calib_status
  
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
    query += " order by enter_time desc where ";
    query += "(completion = 'OK') and instrument ='"; query += instrument;
    query += "' and calib_type ='";
    query += calibType;
    query += "' and '";
    query += timestamp.timeString();
    query += "'> vstart) and (vend > '";
    query += timestamp.timeString();
    query += "')";
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

      int ret = mysql_query(readCxt, q.c_str());
      if (ret) {
        std::cerr << "MySQL error during SELECT, code " << ret << std::endl;
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

  bool Metadata::getReadInfo(unsigned int serialNo, 
                             std::string* dataFmt, 
                             std::string* fmtVersion,
                             std::string* filename) {

    std::string q("select data_fmt, fmt_version, data_ident from ");
    q += table;
    q += "where ser_no='";
    std::strstream s;
    s << serialNo;
    q += s.str(); q+= "'";

    int ret = mysql_query(readCxt, q.c_str());
    if (ret) {
      std::cerr << "MySQL error during SELECT, code " << ret << std::endl;
      return RETMySQLError;
    }
      
    MYSQL_RES *myres = mysql_store_result(readCxt);
    
    if (mysql_num_rows(myres) ) {  // must have been a good serial number
      MYSQL_ROW myRow = mysql_fetch_row(myres);
      *dataFmt = std::string(myRow[0]);
      *fmtVersion = std::string(myRow[1]);
      *filename = std::string(myRow[2]);
      return true;
    }
    else return false;
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
    else if ((testBit = LEVELSuper) & *levelMask) {
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
  bool Metadata::checkCompletionInput(const std::string& stat) {
    return ((stat.compare("OK") == 0 ) || (stat.compare("INCOMPLETE") == 0) 
            || (stat.compare("ABORTED") == 0) );
  }
  bool Metadata::checkProcLevelInput(const std::string& level) {
    return ((level.compare("TEST") == 0) || 
            (level.compare("DEVELOPMENT") == 0) ||
            (level.compare("PRODUCTION") == 0) || 
            (level.compare("SUPERSEDED") == 0 )    );
  }
}
