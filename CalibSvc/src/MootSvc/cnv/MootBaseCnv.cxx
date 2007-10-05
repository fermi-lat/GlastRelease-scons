// $Header$

#include "MootBaseCnv.h"

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
// #include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
// #include "GaudiKernel/IConversionSvc.h"
//#include "GaudiKernel/MsgStream.h"
#include "../MootSvc.h"
#include "CalibSvc/MootTds.h"

MootBaseCnv::MootBaseCnv( ISvcLocator* svc, const CLID& clid) :
  Converter(MOOT_StorageType, clid, svc),
  m_mootSvc (0), m_q(0), m_log(0) {}

StatusCode MootBaseCnv::initialize() {
  StatusCode sc = Converter::initialize();
  if (sc != StatusCode::SUCCESS) return sc;



  // Why am I even doing this?  
  /*
  serviceLocator()->getService ("CalibDataSvc",
                                IID_IDataProviderSvc,
                                (IInterface*&)dp);
  */
  // get mootSvc
  serviceLocator()->getService("MootSvc",  IID_IMootSvc, 
                               (IInterface*&)m_mootSvc);

  m_q = m_mootSvc->getConnection();
  if (m_q) return StatusCode::SUCCESS;
  else return StatusCode::FAILURE;
}

StatusCode MootBaseCnv::finalize() {
  return StatusCode::SUCCESS;   // until we figure out what else to do
}

const unsigned char MootBaseCnv::storageType() {
  return MOOT_StorageType;
}

// Normally this should not get called;  Nodes are not of this type,
// only of derived types
StatusCode MootBaseCnv::createObj(IOpaqueAddress* addr,
                                  DataObject*& refpObject) {
  // From addr get information for type and subtype
  MOOTTYPE tp = (MOOTTYPE) addr->ipar()[0];
  MOOTSUBTYPE sub = (MOOTSUBTYPE) addr->ipar()[1];
  refpObject = new MootBase(tp, sub);
  return StatusCode::SUCCESS;
}


