// $Header$
#ifndef CALIBUTIL_METADATA_H
#define CALIBUTIL_METADATA_H

//#include "calibUtil/Timestamp.h"
#include "facilities/Timestamp.h"
//#include <my_global.h>
//#include <mysql.h>

#ifdef  WIN32
#include <windows.h>
#endif

typedef struct st_mysql MYSQL;

namespace calibUtil {
  /** Provide interface between calibration clients and the
   MySQL database for calibration metadata.  Supported operations
   include writing a new record, looking for the serial number of
   a "best match" record (findBest), and a method for retrieving 
   the information necessary to read the data set corresponding to
   a particular metadata record (getReadInfo). Writing a record 
   is done in stages:  open the record (openRecord), add certain 
   fields (addValidInterval, addNotes, etc.), and finally to write
   the full record to the database (insertRecord).

   The class contains definitions for several enumerated types, 
   corresponding to possible entries for columns in a metadata record.
   Within the MySQL database the values are kept as strings.  The
   enumerated types are provided to insulate clients from these 
   details; enumerations are easier to document and to check for validity.
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

    enum eInstrument {
      INSTLat = 0,
      INSTBtem,
      INSTBfem,
      INSTEm,
      INSTFu,
      INSTCu
    };

    /**
       Calibration types which might conceivably be maintained by this
       system.  We might not actually need them all, but numbers are cheap.
       Note some of these types correspond to more than one entry in
       various tables of calibration types, e.g. the one in
   http://www-glast.slac.stanford.edu/LAT/INT/SVAC/Calibration/calibrations.htm
       For example CTYPE_TKRBadChan handles both hot and dead (but any
       individual dataset has only all hot or all dead).  Similarly for
       CTYPE_TKRAlign, which may apply to any of several subcategories of
       alignment data, and CTYPE_ACDThresh (veto threshold or high threshold).
    */

    enum eCalibType {
      CTYPE_ACDEff = 0,
      CTYPE_ACDThresh,
      CTYPE_ACDPed,
      CTYPE_ACDElecGain,
      CTYPE_TKRBadChan = 8,
      CTYPE_TKRThresh,
      CTYPE_TKRTotSig,
      CTYPE_TKRTotCount,
      CTYPE_TKRAlign,
      CTYPE_CALLightAtt = 16,
      CTYPE_CALLightAsym,
      CTYPE_CALLightYield,
      CTYPE_CALScintEff,
      CTYPE_CALPed,
      CTYPE_CALElecGain,
      CTYPE_CALIntNonlin,
      CTYPE_CALDiffNonlin
    };
      
    /// Constructor keeps track of table of interest
    Metadata(std::string host="*",
             std::string table="*");

    ~Metadata();


    /// Add setting of creator column to row-in-progress
    eRet addCreator(const std::string& creator);

    /// Add description of input to cal procedure to row-in-progress
    eRet addInputDesc(const std::string& desc);

    /// Add notes column to row-in-progress
    eRet addNotes(const std::string& notes);

    /// Set validity interval: period over which calibration data
    /// is applicable.
    eRet addValidInterval(facilities::Timestamp startTime, 
                          facilities::Timestamp endTime);

    /// Add flavor column to row-in-progress.  If this method is not
    /// invoked, the new row will have the default value for flavor,
    /// "VANILLA"
    eRet addFlavor(const std::string& flavor);

    /** Explicit clear of record.  Must have a call to either insertRecord
     *  (to actually write the record to the database) or clearRecord 
     *  (to abort) between successive calls to openRecord.
     */
    void clearRecord();


    /** Return serial number for calibration which is best match to
        criteria, expressed using enumerated types for calibType and
        instrument.  Most clients should use this version of findBest.
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
                 @param flavor       optionally specify non-standard
                                     calibration flavor 
                 @return             status. Should be RETOk.
                                     

       If there are multiple calibrations which are not distinguished
       by the above, pick the one most recently written.
    */
    eRet findBest(unsigned int *ser,
                  eCalibType calibType, 
                  const facilities::Timestamp& timestamp,
                  unsigned int levelMask, 
                  eInstrument instrument,
                  const std::string& flavor="VANILLA");
    /** Return serial number for calibration which is best match to
        criteria, using strings for calibType and instrument arguments.
        This method may be useful for development when a particular
        instrument or calibration type is not officially supported.
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
                 @param flavor       optionally specify non-standard
                                     calibration flavor 
                 @return             status. Should be RETOk.
                                     

       If there are multiple calibrations which are not distinguished
       by the above, pick the one most recently written.
    */
    eRet findBest(unsigned int *ser,
                  const std::string& calibType, 
                  const facilities::Timestamp& timestamp,
                  unsigned int levelMask, 
                  const std::string& instrument,
                  const std::string& flavor="VANILLA");

    const std::string* const getCalibTypeStr(eCalibType cType);
    const std::string* const getDataFmtStr(eDataFmt fmt);
    const std::string* const getInstrumentStr(eInstrument inst);


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
                        
    /** Write a record to the metadata database. Any required columns
     *  not specified by caller will be set to default values.
     @param serialNo   serial number of new record if insert succeeds [output]
     @return           See eRet enumerated type for possible values
     */
    eRet insertRecord(unsigned int* serialNo = 0);

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
    eRet openRecord(eInstrument  inst,
                    eCalibType   ctype,
                    eDataFmt     fmt,
                    const std::string& fmtVersion,
                    const std::string& dataIdent, 
                    eCompletion  completion,
                    eLevel       procLevel = LEVELTest);

    eRet openRecord(const std::string& instr, 
                    const std::string& calibType,
                    const std::string& dataFmt,
                    const std::string& fmtVersion,
                    const std::string& dataIdent, 
                    const std::string& completion,
                    const std::string& level);



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

    MYSQL* m_readCxt;
    MYSQL* m_writeCxt;

    bool addLevel(std::string& q, unsigned int *levelMask);

    static bool connect(MYSQL * cxt, const std::string& user, 
                        const std::string& pw, eRet& err);

    bool connectRead(eRet& err);
    bool connectWrite(eRet& err);

    //    static void makeQuery(std::string& query, unsigned int *levelMask);
    const std::string*  const getCompletionStr(eCompletion cmp);
    const std::string* const getProcLevelStr(eLevel level);

    

    /// Discover username and add to row-in-progress
    eRet addUser();

    /** Keep track of which columns in row have been initialized 
        with bit mask */
    enum eRow {
      eOpened = 1,
      eValid = 2,
      eInputDesc = 4,
      eComment = 8,
      eFlavor  = 0x10,
      eCreator = 0x20 };

    /// Constant bit mask indicating all necessary fields have been set
    static const unsigned int s_rowReady;

    std::string m_row;     // place to keep row as it's being built
    unsigned int m_rowStatus;
    std::string  m_host;
    std::string  m_table;
  };
}

#endif
