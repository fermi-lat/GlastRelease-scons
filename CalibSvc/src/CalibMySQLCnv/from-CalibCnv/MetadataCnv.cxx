/**
  @file MetadataCnv.cxx

  @author Joanne Bogart
  $Header$
*/

#include "MetadataCnv.h"
#include "CalibData/MetadataInfo.h"

// Following not needed since MetadataCnv.h references Converter.h
// which references IConverter.h which references ClassID.h
// #include "GaudiKernal/ClassID.h"   // for storage type id's

MetadataCnv::MetadataCnv(ISvcLocator* svc) : m_meta(0) {
  using CalibData::MetadataInfo;
  //  Converter(MYSQL_StorageType, CLID_Calib_MetadataInfo, svc);
  Converter(MYSQL_StorageType, MetadataInfo::classID(), svc);
  // can't think of anything else to do here
}

StatusCode MetadataCnv::initialize() {
  StatusCode status = Converter::initialize();

  // Locates the Detector Data Service 
  IDataProviderSvc* dp;
  serviceLocator()->getService ("DetectorDataSvc",
                                IID_IDataProviderSvc,
                                (IInterface*&)dp);
  setDataProvider(dp);
  
  // Locate the calib Conversion Service
  serviceLocator()->getService ("CalibCnvSvc",
                                IID_ICalibCnvSvc,
                                (IInterface*&)m_CalibCnvSvc);

  // Locate det data service
  serviceLocator()->getService("DetectorDataSvc",
                               IID_IDetDataSvc,
                               (IIInterface*&)m_DetDataSvc);

  // returns
  return status;
}

StatusCode MetadataCnv::finalize() {
  // RIP dear grand father!
  return Converter::finalize();
}


StatusCode MetadataCnv::createObj(IOpaqueAddress* addr,
                                  DataObject*& refpObject) {
   // creates a msg stream for debug purposes
   MsgStream log( msgSvc(), "MetadataCnv" );
   
   // maked sure the address is not null
   if (0 == addr) {
     return StatusCode::FAILURE;
   }



}


StatusCode MetadataCnv::updateObj(IOpaqueAddress* addr,
                                  DataObject*& refpObject) {

}
 
