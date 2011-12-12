// $Header$
#ifndef CalibDataSvc_h
#define CalibDataSvc_h

// Base classes
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/Time.h"

#include "CalibSvc/IInstrumentName.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/CalibTime.h"
// Forward declarations
//class ITime;
class StatusCode;
class IDataProviderSvc;
class MsgStream;

/** @class CalibDataSvc

    A DataSvc specialized for calibration data.  
    This Service borrows heavily from DetDataSvc.  In particular it
    implements the IDetDataSvc interface.  The only significant
    difference is in initialize() and in the elimination of members
    concerned with detector (geometry) description.

    Maybe will also need to implement another abstract service
    which gets (and sets?) instrument.

    @author J. Bogart
    @date   15 Oct. 2002   
*/
class IAddressCreator;
class MsgStream;

class CalibDataSvc  : public DataSvc,
                      virtual public IDetDataSvc,
                      virtual public IIncidentListener,
                      virtual public IInstrumentName,
                      virtual public ICalibPathSvc
{    

  friend class SvcFactory<CalibDataSvc>;

public:
  
  // Overloaded DataSvc methods

  virtual StatusCode initialize();
  
  virtual StatusCode finalize();
  
  /// Remove all data objects in the data store.
  virtual StatusCode clearStore();
  
  /// Update object
  virtual StatusCode updateObject( DataObject* toUpdate );

  /// Load object.  Override DataSvc implementation to get current 
  /// event time first if necessary
  virtual StatusCode loadObject(IConversionSvc* pLoader, IRegistry* pRegistry);


protected:
  
  CalibDataSvc(const std::string& name, ISvcLocator* svc);
  virtual ~CalibDataSvc();
  
public:
  
  // Reimplemented from IInterface

  /// Query the interface of the service
  virtual StatusCode queryInterface( const InterfaceID& riid, 
				     void** ppvInterface );  

public:

  // Implementation of the IDetDataSvc interface

  /// Check if the event time has been set 
  virtual bool validEventTime() const ;

  /// Get the event time  
  virtual const Gaudi::Time& eventTime() const;

  /// Set the new event time  
  virtual void setEventTime(const Gaudi::Time& time);

public:
  //Implementation of IInstrumentName interface

  /// Check if the event time has been set 
  virtual const bool validInstrumentName() const;

  /// Get the instrument name
  virtual const std::string& getInstrumentName() const;

  /// Set the instrument name
  virtual void setInstrumentName(const std::string& name);

  // Implementation of ICalibPathSvc interface
public:

  virtual const std::string
  getCalibPath(const ICalibPathSvc::CalibItem item, 
               const std::string& flavor="") const;

  
  // Implementation of the IIncidentListener interface

  /// Inform that a new incident has occured
  virtual void handle( const Incident& );

  /// For use of CalibMySQLCnvSvc, to set "use event time mode"
  virtual void setUseEventTime(bool useEventTime) {
    m_useEventTime = useEventTime;}

 private:
  //properties
  /// Calibration Data Persistency Storage type 
  long int              m_calibStorageType;


  /// Name of the root node of the calib store
  std::string      m_calibRootName;

  /// calibration types
  StringArrayProperty m_calibList;

  //  StringArrayProperty m_flavorList;
  std::vector<std::string>  m_flavorList;

  /// Has the event time been defined?
  bool             m_eventTimeDefined;

  /// Current event time
  Gaudi::Time*           m_eventTime; 

  bool             m_instrumentDefined;
  std::string      m_instrumentName;
  std::string      m_timeSource;

  /// Just for diagnostic purposes, keep count of events seen
  unsigned m_nEvent;

  /// Set to true by handle( ) at BeginEvent; cleared after timestamp acquired
  bool m_newEvent;

  /// Dumping place for various time-fetching methods to save the timestamp
  CalibData::CalibTime m_time;

  IDataProviderSvc* m_eventSvc;
  /** Redundant job option to indicate whether or not to check for valid
      event time when fetching calibrations.  There already is a similar thing
      for CalibMySQLCnvSvc, but no easy way for one service to tell the 
      other about it.
  */
  bool  m_useEventTime;


  /// Private utility, called from initialize()
  StatusCode makeFlavorNodes(IAddressCreator*  calibCreator);

  /// Private utility to check if internal timestamp has been updated for
  /// the event; if not do it.
  StatusCode updateTime();


  /*
  /// Typedef for ptr to function which will attempt to fetch 
  /// "current time" and store in CalibTime arg.
  typedef StatusCode (*FetchTime)(CalibData::CalibTime& ourTime,
                                  MsgStream& log);

  FetchTime m_fetcher;
  */

  enum TIMESOURCE {
    TIMESOURCEnone = 0,
    TIMESOURCEdata,
    TIMESOURCEmc,
    TIMESOURCEclock,
    TIMESOURCEdigi
  };

  TIMESOURCE m_timeSourceEnum;

  /// Fetch time from mc data (no digis)
  StatusCode fetchMcTime();


  // Supersedes above.  For data, digi and mc *should* always fetch time
  // from the same place.  
  // TODO - determine if I can get rid of fetchMcTime 
  /// Fetch time from Event
  StatusCode fetchEventTime();

  /// Fetch time from fake clock, using parameters below
  StatusCode fetchFakeClockTime();

  /// Initialize path arrays
  void initPathArrays();

  /// Absolute time of first event (yyyy-mm-dd_hh:mm, trailing fields
  /// optional)
  std::string m_startTimeAsc;

  /// Absolute time of first event (seconds)
  long m_startTime;

  /// Absolute time spacing between events (milliseconds)
  long m_delayTime;

  unsigned m_LATCMaster;
  MsgStream *m_log;

  /// Make it easy to retrieve path (not including flavor) from enum
  /// Might later include config paths here, or they might have their
  /// own vector
  std::vector<CLID>    m_calibCLIDs;
  std::vector<std::string> m_calibPaths;
};

#endif //  CalibDataSvc_h



