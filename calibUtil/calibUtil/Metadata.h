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
      RETBadValue = 5,
      RETMySQLError = 6
    };
    /// Used to form bit masks for dbs queries
    enum eLevel {
      LEVELProd = 1,
      LEVELDev  = 2,
      LEVELTest = 4,
      LEVELSuperseded = 8
    };

    enum eDataFmt {
      FMTXml = 0,
      FMTRoot = 1,
      FMTUnknown = 2
    };

    enum eCompletion {
      CMPLOk = 0,
      CMPLInc = 1,
      CMPLAbort = 2,
      CMPLUnknown = 3
    };
      
    /// Constructor keeps track of table of interest
    Metadata();

    ~Metadata();


    /** Start a new metadata record by supplying all absolutely
        required information as arguments:
        @param instr  Instrument name, such as "flight", "EM",etc.
        @param calibType   One of the recognized calibration types,
                           such as "hotStrips"
        @param dataFmt  For now, one of "XML" or "ROOT"
        @param fmtVersion  Something to further identify data format
                           in case it evolves with time, e.g., "f1v0"
        @param dataIdent   Identifies the data being described in the
                           record; typically the full file spec of 
                           the file containing the data
        @param completion  Completion status of calibration; has nothing
                           to do with health of the detector being 
                           calibrated.  Possible values are "OK", 
                           "INC" or "ABORT"
        @return            See the eRet enumerated type for possible
                           values
    */
    eRet openRecord(const std::string& instr, 
                    const std::string& calibType,
                    eDataFmt     fmt,
                    const std::string& fmtVersion,
                    const std::string& dataIdent, 
                    eCompletion  completion,
                    eLevel       procLevel = LEVELTest);

    /** Write a record to the metadata database. Any required columns
     *  not specified by caller will be set to default values.
     */
    eRet insertRecord();

    /** Explicit clear of record.  Must have a call to either insertRecord
     *  (to actually write the record to the database) or clearRecord 
     *  (to abort) between successive calls to openRecord.
     */
    void clearRecord();

    /// Set validity interval: period over which calibration data
    /// is applicable.
    eRet addValidInterval(Timestamp startTime, Timestamp endTime);

    /// Add setting of creator column to row-in-progress
    eRet addCreator(std::string creator);

    /// Add notes column to row-in-progress
    eRet addNotes(std::string notes);

    /// Add description of input to cal procedure to row-in-progress
    eRet addInputDesc(std::string desc);


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
    eRet findBest(unsigned int *ser,
                  const std::string& calibType, 
                  Timestamp timestamp,
                  unsigned int levelMask, 
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
    eRet getReadInfo(unsigned int serialNo, 
                     eDataFmt&     dataFmt,
                     std::string& fmtVersion,
                     std::string& dataIdent);
                        
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
    void disconnectRead();
    void disconnectWrite();

  private:
    // may change MYSQL structures to not be static
    MYSQL* readCxt;
    MYSQL* writeCxt;

    bool connectRead(eRet& err);
    bool connectWrite(eRet& err);
    static bool connect(MYSQL * cxt, const std::string& user, 
                        const std::string& pw, eRet& err);

    //    static void makeQuery(std::string& query, unsigned int *levelMask);
    bool addLevel(std::string& q, unsigned int *levelMask);
    const std::string*  const checkCompletionInput(eCompletion cmp);
    const std::string* const checkProcLevelInput(eLevel level);
    const std::string* const checkDataFmtInput(eDataFmt fmt);

    

    /// Discover username and add to row-in-progress
    eRet addUser();

    /** Keep track of which columns in row have been initialized 
        with bit mask */
    enum eRow {
      eOpened = 1,
      eValid = 2,
      eInputDesc = 4,
      eComment = 8,
      eCreator = 0x10 };

    /// Constant bit mask indicating all necessary fields have been set
    static const unsigned int rowReady;

    std::string row;     // place to keep row as it's being built
    unsigned int rowStatus;
    std::string  table;
  };
}

#endif
