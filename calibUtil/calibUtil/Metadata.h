// $Header$
#ifndef CALIBUTIL_METADATA_H
#define CALIBUTIL_METADATA_H

#include <my_global.h>
#include <mysql.h>
#include <string>

namespace calibUtil {
  /** Provide interface between calibration clients and the
   MySQL database for calibration metadata.
  */
  class Metadata {
  public:
    enum eError {
      
    };
    /// Used to form bit masks for dbs queries
    enum eDisposition {
      DISPprod = 1,
      DISPdev  = 2,
      DISPtest = 4,
      DISPsuper = 8
    };
      

    /// Constructor initializes connection pointers to "not connected"
    /// value.
    Metadata();

    ~Metadata();


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
    static bool openRecord(const std::string& instr, 
                           const std::string& calibType,
                           const std::string& dataFormat, 
                           const std::string& fmtVersion,
                           const std::string& filename, 
                           const std::string& calibStatus);
    static bool writeRecord(eError& err);
    static void clearRecord();
    static void addField(cosnt std::string& fieldName, 
                         const std::string& fieldValue);

    /** Return serial number for calibration which is best match to
        criteria
                 @param calibType    type of data, must match
                 @param dispMask     acceptable dispositions ("production"
                                     better than "dev" better than "test"
                                     better than "superseded")
                 @param timestamp    must be within validity interval; 
                                     closer to center is better
                 @return serial number of best match as integer or 
                                       zero if no matches.

       If there are multiple calibrations which are not distinguished
       by the above, pick the one most recently written.
    */
    static unsigned int findBest(const std::string& calibType, 
                                 unsigned dispMask, 
                                 unsigned int timestamp);

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
                        

    // Might make these private
    static bool connectRead(eError& err);
    static bool connectWrite(eError& err);
    static bool disconnectRead();
    static bool disconnectWrite();

  private:
    static MYSQL* readCxt;
    static MYSQL* writeCxt;
    static std::string row;     // place to keep row as it's being built
  };
}

#endif
