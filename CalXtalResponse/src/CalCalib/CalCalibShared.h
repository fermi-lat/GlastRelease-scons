#ifndef CalCalibShared_h
#define CalCalibShared_h

// LOCAL INCLUDES
#include "IdealCalCalib.h"

// GLAST INCLUDES
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"


// EXTLIB INCLUDES
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IDataProviderSvc.h"

// STD INCLUDES

/** \brief store items shared amongst several classes which compose CalCalibSvc


*/
class CalCalibShared {
 public:
  CalCalibShared() :
    m_dataProviderSvc(0),
    m_service(0) {}

  /// post construction: run during Gaudi init phase.
  StatusCode initialize(Service &service);

  /// stores 'ideal' flavor generic calibration constants
  IdealCalCalib m_idealCalib;

  IDataProviderSvc *m_dataProviderSvc;

  Service *m_service;

  ///-- jobOptions --///

  /// xml file contains 'ideal' flavor parameters
  StringProperty m_idealCalibXMLPath;

  /// name of CalibDataSvc, main data source
  StringProperty m_calibDataSvcName;     


};

#endif
