// $Header$
#ifndef CalibDataSvc_h
#define CalibDataSvc_h

// Base classes
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/IIncidentListener.h"

#include "CalibSvc/IInstrumentName.h"

// Forward declarations
class ITime;
class StatusCode;

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

class CalibDataSvc  : public DataSvc,
                      virtual public IDetDataSvc,
                      virtual public IIncidentListener,
                      virtual public IInstrumentName
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

protected:
  
  CalibDataSvc(const std::string& name, ISvcLocator* svc);
  virtual ~CalibDataSvc();
  
public:
  
  // Reimplemented from IInterface

  /// Query the interface of the service
  virtual StatusCode queryInterface( const IID& riid, 
				     void** ppvInterface );  

public:

  // Implementation of the IDetDataSvc interface

  /// Check if the event time has been set 
  virtual const bool validEventTime() const ;

  /// Get the event time  
  virtual const ITime& eventTime() const ;

  /// Set the new event time  
  virtual void setEventTime(const ITime& time);

public:
  //Implementation of IInstrumentName interface

  /// Check if the event time has been set 
  virtual const bool validInstrumentName() const;

  /// Get the instrument name
  virtual const std::string& getInstrumentName() const;

  /// Set the instrument name
  virtual void setInstrumentName(const std::string& name);

public:
  
  // Implementation of the IIncidentListener interface

  /// Inform that a new incident has occured
  virtual void handle( const Incident& );

 private:
  //properties
  /// Calibration Data Persistency Storage type
  int              m_calibStorageType;

  /// Name of the root node of the calib store
  std::string      m_calibRootName;

  /// calibration types
  StringArrayProperty m_calibList;
  /*
  /// Flag to control if the persistency is required
  bool             m_usePersistency;
  */

  /// Has the event time been defined?
  bool             m_eventTimeDefined;

  /// Current event time
  ITime*           m_eventTime; 

  bool             m_instrumentDefined;
  std::string      m_instrumentName;

};

#endif //  CalibDataSvc_h



