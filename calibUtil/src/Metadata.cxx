// $Header$


#include "calibUtil/Metadata.h"
#include "facilities/Util.h"
#include "mysql.h"

namespace calibUtil {
  MYSQL* Metadata::readCxt = 0;
  MYSQL* Metadata::writeCxt = 0;
  std::string row("");

  bool Metadata::connectRead(eError& err) {

    return connect(readCxt, std::string("reader"), err);

  }

  bool Metadata::connectWrite(eError& err) {

    return connect(writeCxt, std::string("calibrator"), err);

  }


  bool Metadata::connect(MYSQL* cxt, std::string grp, eError& err)  {
    std::string cnfFile("$(CALIBUTILROOT)/src/.mysqlopt.cnf");

    // Translate environment variables
    int nSub = facilities::Util::expandEnvVar(&cnfFile);

    if (nSub != 1) {
      err = ERRbadCnfFile;
      return false;
    }

    std::string host("$(MYSQL_HOST)");

    nSub = facilities::Util::expandEnvVar(&host);
    if (!nSub) {
      err = ERRbadHost;
      return false;
    }


    mysql_init(cxt);

    mysql_options(cxt, MYSQL_READ_DEFAULT_FILE, cnfFile.c_str());
    mysql_options(cxt, MYSQL_READ_DEFAULT_GROUP, grp.c_str());

    MYSQL *connected = mysql_real_connect(cxt, host.c_str(), 0, 0, "calib",
                                 0, (char *) NULL, 0);
    if (connected != 0) {
      cxt = connected;
      err = ERRok;
      return true;
    }
    else {
      err = ERRnoConnect;
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

  bool Metadata::openRecord(const std::string& instr, 
                            const std::string& calibType,
                            const std::string& dataFormat, 
                            const std::string& fmtVersion,
                            const std::string& filename, 
                            const std::string& calibStatus){
    if (row.size() > 0) {
      std::cout << 
        "calibration metadata record build already in progress" << std::endl;
      return false;
    }

    //    += is a synonym for row.append(..)
    row += "insert into metadata set instrument=";    row += instr;
    row += ", calib_type=";  row += calibType;
    row += ", data_format="; row += dataFormat;
    row += ", format_version="; row += fmtVersion;
    row += ", calib_status="; row += calibStatus;
    row += ", calib_file_name"; row += filename;
    return true;
  }

  /*

| Field               | Type           | Null | Key | Default | Extra         |
+---------------------+----------------+------+-----+---------+---------------+
| serial_num          | int(11)        |      | PRI | NULL    | auto_increment|
| instrument          | varchar(20)    | YES  |     | NULL    |  [open]       |
| calib_type          | varchar(40)    | YES  |     | NULL    |  [open]       |
| calib_type_num      | int(11)        | YES  |     | NULL    |               |
| data_format         | varchar(20)    |      |     |         |  [open]       |
| format_version      | varchar(20)    | YES  |     | NULL    |  [open]       |
| data_id             | varchar(40)    | YES  |     | NULL    |               |
| data_size           | int(11)        | YES  |     | NULL    |               |
| validity_start_time | datetime       | YES  |     | NULL    |               |
| validity_end_time   | datetime       | YES  |     | NULL    |               |
| proc_time           | datetime       | YES  |     | NULL    |               |
| proc_level          | enum('TEST','DEVELOPMENT','PRODUCTION','SUPERSEDED')
                                       | YES  |     | NULL    |               |
| calib_status        | enum('OK','INCOMPLETE','ABORTED')        [open]
                                       | YES  |     | NULL    |               |
| creator             | varchar(40)    | YES  |     | NULL    |               |
| input_desc          | varchar(255)   | YES  |     | NULL    |               |
| calib_file_name     | varchar(255)   |      |     |         |  [open]       |
| comments            | varchar(255)   | YES  |     | NULL    |               |
+---------------------+----------------+------+-----+---------+---------------+

  */


  bool Metadata::insertRecord(eError& err) {
    // To be written
    return true;
  }

  // Explicit clear of in-progress record building
  void Metadata::clearRecord() {
    row.resize(0);
  }

  void Metadata::addField(const std::string& fieldName, 
                          const std::string& fieldValue) {
    // To be written

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
