// $Header$


#include "calibUtil/Metadata.h"
#include "facilities/Util.h"
#include "mysql.h"

namespace calibUtil {
  MYSQL* Metadata::readCxt = 0;
  MYSQL* Metadata::writeCxt = 0;
  std::string Metadata::row("");
  const unsigned int Metadata::rowReady = Metadata::eOpened | Metadata::eValid 
  | Metadata::eInputDesc | Metadata::eComment;

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

  Metadata::eRet Metadata::openRecord(const std::string& instr, 
                                      const std::string& calibType,
                                      const std::string& dataFormat, 
                                      const std::string& fmtVersion,
                                      const std::string& filename, 
                                      const std::string& calibStatus,
                                      const std::string& procLevel){
    if (row.size() > 0) {
      std::cout << 
        "calibration metadata record build already in progress" << std::endl;
      return RETWrongState;
    }

    //    += is a synonym for row.append(..)
    row += "insert into metadata set instrument=";    row += instr;
    row += ", calib_type=";  row += calibType;
    row += ", data_format="; row += dataFormat;
    row += ", format_version="; row += fmtVersion;
    row += ", calib_status="; row += calibStatus;
    row += ", calib_file_name"; row += filename;
    row += ", proc_level"; row += procLevel;
    rowStatus = eOpened;
    return RETOk;
  }

  /*

| Field               | Type           | Null | Key | Default | Extra         |
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


  Metadata::eRet Metadata::insertRecord() {
    // Send it off (To be written)

    if ((!(rowStatus & rowReady)) == rowReady) return RETWrongState;

    if (!(rowStatus & eCreator)) {
      addCreator("calibUtil::Metadata::insertRecord");
    }
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
    
    row += ", validity_start_time="; row += startTime.timeString();
    row += ", validity_end_time="; row += endTime.timeString();
    rowStatus |=  eValid;
    return RETOk;   // or something else
  }

  Metadata::eRet Metadata::addComment(std::string comment) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eComment) return RETWrongState;

    row += ",comments="; row += comment;
    rowStatus |= eComment;
    return RETOk;
  }

  Metadata::eRet Metadata::addInputDesc(std::string desc) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eInputDesc) return RETWrongState;

    row += ",input_desc="; row += desc;
    rowStatus |= eInputDesc;
    return RETOk;
  }

  Metadata::eRet Metadata::addCreator(std::string creator) {
    if (!(rowStatus & eOpened) ) return RETWrongState;
    if (rowStatus & eCreator) return RETWrongState;

    row += ",creator="; row += creator;
    rowStatus |= eCreator;
    return RETOk;
  }

  unsigned int Metadata::findBest(const std::string& calibType, 
                                  unsigned dispMask, 
                                  unsigned int timestamp) {
    // To be written
    return 0;

  }

  bool Metadata::getReadInfo(unsigned int serialNo, 
                             const std::string* dataFormat, 
                             const std::string* fmtVersion,
                             const std::string* filename) {
    // To be written
    return true;
  }
           
}
