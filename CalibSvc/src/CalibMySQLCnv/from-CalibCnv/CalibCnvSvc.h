// $Header$
#ifndef CALIBCNV_CALIBCNVSVC_H
#define CALIBCNV_CALIBCNVSVC_H

#include "CalibCnvSvc/ICalibSvc.h"
#include "GaudiKernel/ConversionSvc.h"
// #include "GaudiKernel/ICnvManager.h"
// Probably need the following:
#include "GaudiKernel/GenericAddress.h"     

// Forward and external declarations
template <class TYPE> class SvcFactory;

namespace calibUtil {
  class Metadata;
}

// Don't think I need these
// class IConversionSvc;
// class IDetDataSvc;

/** @class CalibCnvSvc 
 *
 *  A conversion service for calibration data. It is based on the generic
 *  conversion service and implements in addition the ICalibSvc interface.
 *
 *  @author Joanne Bogart
 */
class CalibCnvSvc : public ConversionSvc,
                    virtual public ICalibSvc {

  /// Friendship needed to access protected constructor
  friend class SvcFactory<CalibCnvSvc>;
  
protected:
  CalibCnvSvc(const std::string& name, ISvcLocator* svc);

  virtual ~CalibCnvSvc();

public:
  // From IInterface
  
  /**
     Queries interfaces of Interface.
     @param riid ID of Interface to be retrieved
     @param ppvInterface Pointer to Location for interface pointer
     @return status depending on the completion of the call
   */
  virtual StatusCode queryInterface (const IID& riid, void** ppvInterface);
  
  // Overload ConversionSvc (from Service) initialize
  virtual StatusCode initialize();

  // Overload ConversionSvc (from Service) finalize
  virtual StatusCode finalize();

  // Most likely can get by with ConversionSvc implementation of
  // createObj, fillObjRefs, updateObj, updateObjRefs 
  // It just calls the correct converter for the particular
  // object for each of these functions.  
  // Or, might want to include empty implemntation for fillObjRefs and 
  // updateObjRefs since we don't expect to have any refs in our
  // objects.


  /**
   * Create a Calib address using explicit arguments to identify a 
   * single object
   * @param svc_type the service type
   * @param CLID the CLID of the calib. data for which an address is created
   * @param par an array of two strings containing the instrument name
   * and objectname in this order
   * @param ip can be ignored 
   * @param refpAddress the new address created
   * @return a StatusCode giving the status of the address creation
   */
  virtual StatusCode createAddress(unsigned char svcType,
                                   const CLID& clid,
                                   const std::string* stringPar, 
                                   const unsigned long* longPar,
                                   IOpaqueAddress*& refpAddress) {

  }

  // Implementation of ICalibCnvSvc
  virtual StatusCode findBest(std::string calibType,
                              std::string instrument,
                              facilities::Timestamp* eventTime,
                              CalibKey *key);

  virtual StatusCode locatePDS(CalibKey key, std::string& dataFmt,
                               std::string& fmtVersion,
                               std::string& dataIdent);
                   

protected:
  /**
   * Use a factory so that standard constructor need not be exposed
   * @param name String with service name
   * @param svc Pointer to service locator interface
   * @return Reference to CalibCnvSvc
   */
  CalibCnvSvc(const std::string& name, ISvcLocator* svc);
  
  virtual ~CalibCnvSvc();

private:

  calibUtil::Metadata*  m_meta;
  // Bunch of stuff will go here ultimately

  /*   Here is the stuff we probably don't need from
     LHCb-cond/Det/DetCond/src/component/ConditionsDbCnvSvc.h

  /// Global tag name (can be set using the JobOptionsSvc)
  std::string          m_globalTag;

  /// Handle to the low-level gate to the ConditionsDB
  ConditionsDBGate*    m_conditionsDBGate;

  /// Handle to the IConversionSvc interface of the DetectorPersistencySvc
  IConversionSvc*      m_detPersSvc;

  
  */
  /// Handle to the IDetDataSvc interface of the (calib) data service
  ICalibDataSvc*         m_detDataSvc;


};
#endif
