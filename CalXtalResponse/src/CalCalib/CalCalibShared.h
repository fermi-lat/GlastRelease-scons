#ifndef CalCalibShared_h
#define CalCalibShared_h

// LOCAL INCLUDES
#include "IdealCalCalib.h"

// GLAST INCLUDES
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"


// EXTLIB INCLUDES
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "CalibSvc/ICalibPathSvc.h"

// STD INCLUDES

/** \brief store items shared amongst several classes which compose CalCalibSvc,
    \note it would be overkill to make all the implementation sub-classes Gaudi classes


*/
class CalCalibShared {
 public:
  CalCalibShared() :
    m_dataProviderSvc(0),
    m_calibPathSvc(0),
    m_service(0) {}

  /// post construction: run during Gaudi init phase.
  StatusCode initialize(Service &service);

  /// xml file contains 'ideal' flavor parameters
  StringProperty m_idealCalibXMLPath;

  /// stores 'ideal' flavor generic calibration constants
  IdealCalCalib m_idealCalib;

  /// name of CalibDataSvc, main data source
  StringProperty m_calibDataSvcName;     

  /// CalibDataSvc
  IDataProviderSvc *m_dataProviderSvc;

  /// CalibPathSvc (should be provided by CalibDataSvc
  ICalibPathSvc*  m_calibPathSvc;

  /// allow me to look up services
  Service *m_service;


};

#endif
