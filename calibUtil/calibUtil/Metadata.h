// $Header$
#ifndef CALIBUTIL_METADATA_H
#define CALIBUTIL_METADATA_H

#include "calibUtil/Timestamp.h"
//#include <my_global.h>
//#include <mysql.h>

typedef struct st_mysql MYSQL;

namespace calibUtil {
  /** Provide interface between calibration clients and the
   MySQL database for calibration metadata.  There is no need for
   anything other than a bunch of static routines and a little
   bit of private static data to keep track of connections, etc.
  */
  class Metadata {
  public:
    enum eRet {
      RETOk = 0,
      RETBadCnfFile = 1,
      RETBadHost = 2,
      RETNoConnect = 3,
      RETWrongState = 4,
      RETMySQLError = 5
    };
    /// Used to form bit masks for dbs queries
    enum eLevel {
      LEVELProd = 1,
      LEVELDev  = 2,
      LEVELTest = 4,
      LEVELSuper = 8
    };
      

    /// Nothing to do in constructor since everything is static.
    Metadata() {};

    ~Metadata() {};


    /** Start a new metadata record by supplying all absolutely
        required information as arguments:
        @param instr  Instrument name, such as "flight", "EM",etc.
        @param calibType   One of the recognized calibration types,
                           such as "hotStrips"
        @param dataFormat  For now, one of "XML" or "ROOT"
        @param fmtVersion  Something to further identify data format
                           in case it evolves with time, e.g., "f1v0"
        @param filename    Full file spec of file containing the data
        @param calibStatus Completion status of calibration; has nothing
                           to do with health of the detector being 
                           calibrated.  Possible values are "OK", 
                           "INCOMPLETE" or "ABORTED"
        @return     false if there is already an open record; else
                    true
        The rest is optional or has a sensible default. Other fields include
           serial_num (index automatically incremented by MySQL)
           calib_type_num (maintained by calibUtil)
           data_size      (optional; to be set by client if used)
           proc_time      (defaults to current time)
           validity_start_time, validity_end_time  (optional; default
           to current time)

    */
    static eRet openRecord(const std::string& instr, 
                             const std::string& calibType,
                             const std::string& dataFormat, 
                             const std::string& fmtVersion,
                             const std::string& filename, 
                             const std::string& calibStatus,
                             const std::string& procLevel = "TEST");

    /** Write a record to the metadata database. Any required columns
     *  not specified by caller will be set to default values.
     */
    static eRet insertRecord();

    /** Explicit clear of record.  Must have a call to either insertRecord
     *  (to actually write the record to the database) or clearRecord 
     *  (to abort) between successive calls to openRecord.
     */
    static void clearRecord();

    //    static void addField(const std::string& fieldName, 
    //                         const std::string& fieldValue);

    /** Routine for setting validity interval.
     */
    static eRet addValidInterval(Timestamp startTime,
                                 Timestamp endTime);

    static eRet addComment(std::string comment);

    static eRet addInputDesc(std::string desc);

    static eRet addCreator(std::string creator);

    /** Return serial number for calibration which is best match to
        criteria
                 @param ser          serial number of best match 
                                     as integer or zero if no matches
                                     (output)
                 @param calibType    type of data, must match
                 @param timestamp    must be within validity interval; 
                                     closer to center is better
                 @param levelMask    acceptable levels ("production"
                                     better than "dev" better than "test"
                                     better than "superseded")
                 @param instrument   e.g. LAT, EM, CU,...
                 @return             status. Should be RETOk.
                                     

       If there are multiple calibrations which are not distinguished
       by the above, pick the one most recently written.
    */
    static eRet findBest(unsigned int *ser,
                         const std::string& calibType, 
                         Timestamp timestamp,
                         unsigned levelMask, 
                         const std::string& instrument);


    // Might also want a "findAll" which would just return a list
    // of serial numbers, and a "getRecord" which would either
    // just return the full row as a string or parse it into 
    // its separate columns

    /** Given a calibration serial number, return information needed for 
        caller to read in the data.  
          @param  serialNo           [input]
          @param  dataFormat
          @param  fmtVersion
          @param  filename
          @return     true if serialNo exists in dbs and "filename" has
                      non-null value; else false.
    */
    static bool getReadInfo(unsigned int serialNo, 
                            const std::string* dataFormat, 
                            const std::string* fmtVersion,
                            const std::string* filename);
                        
  /** 
    // Additional services will probably be needed to
    //   1 change proc_level of a given calibration, e.g. from
    //     PRODUCTION to SUPERSEDED
    //   2 "split" a calibration.  That is, make two new calibration
    //     metadata records, both pointing to the same actual data
    //     file, with validity periods whose union = validity period
    //     of an existing metadata record pointing to that file.  Then
    //     (perhaps) mark the original metadata record as superseded
    //
    //  Why would one want #2?  If some reprocessing were done, but
    //  just on some sub(time)interval, then the old calibration would
    //  no longer be preferred for the full interval; want to be able
    //  to mark it as superseded for the subinterval, still usable
    //  for the remainder.  
  */

    // Might make these private
    static bool connectRead(eRet& err);
    static bool connectWrite(eRet& err);
    static void disconnectRead();
    static void disconnectWrite();

  private:
    static MYSQL* readCxt;
    static MYSQL* writeCxt;
    static std::string row;     // place to keep row as it's being built
    static bool connect(MYSQL* cxt, const std::string& group, eRet& err);

    //    static void makeQuery(std::string& query, unsigned int *levelMask);
    static bool addLevel(std::string& q, unsigned int *levelMask);

    /** Keep track of which columns in row have been initialized 
        with bit mask */
    enum eRow {
      eOpened = 1,
      eValid = 2,
      eInputDesc = 4,
      eComment = 8,
      eCreator = 0x10 };
    static const unsigned int rowReady;
    static unsigned int rowStatus;
  };
}

#endif
